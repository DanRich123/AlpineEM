from deap import base, creator, tools, algorithms
import numpy as np
import random
import os, csv, time, shutil, subprocess

# ============================================================
# USER SETTINGS
# ============================================================

AUTO_CLEANUP = True           # Delete generation folders and slurm scripts afterward (set False if debugging something or if we want to keep that data)
SEED = False                  # can submit an intial numpy file 'seeds.npy' to replace initial random generation
CHECK_INTERVAL = 15           # Seconds between result checks - too long slows it down, too fast and Alpine yells at you
NGEN = 3                      # Total generations - first is random or seeded - labeled generation 0 or generation initialization
MU = 20                       # Survivors of each generation to take to the next
LAMBDA = 20                   # Offspring per generation that are created
CXPB, MUTPB = 0.7, 0.3        # Probabilities of crossover vs mutation occuring

# Files needed by FDTD simulations
NEEDED_FILES = [
    "FDTD_full_integration_bias1.py",
    "FDTD_full_integration_bias2.py",
    "Post_Processor_bias1.py",
    "Post_Processor_bias2.py",
    "metal.dat",
    "clear.dat",
    "test"
]

# ============================================================
# FITNESS FUNCTION
# ============================================================

def run_program(individual, work_dir="."):
    os.makedirs(work_dir, exist_ok=True)

    # Save parameters
    params = [int(x) for x in individual]
    np.savetxt(os.path.join(work_dir, "params.txt"), params)

    def run(cmd):
        subprocess.run(cmd, cwd=work_dir, shell=True, check=True)

    try:
        # Run bias 1 simulation
        run("python FDTD_full_integration_bias1.py")
        run("python Post_Processor_bias1.py")
        data = np.loadtxt(os.path.join(work_dir, "S_parameters.csv"), delimiter=",", skiprows=1).T
        freq, s11_bias1_phase = data[0], np.unwrap(data[2])

        # Run bias 2 simulation
        run("python FDTD_full_integration_bias2.py")
        run("python Post_Processor_bias2.py")
        data = np.loadtxt(os.path.join(work_dir, "S_parameters.csv"), delimiter=",", skiprows=1).T
        s11_bias2_phase = np.unwrap(data[2])

        # Calculate fitness
        phase_resolution = 15 * np.pi / 180
        phase_diff = s11_bias1_phase - s11_bias2_phase

        count = 0
        for f, p in zip(freq, phase_diff):
            if 5 <= f <= 8 and p < phase_resolution:
                count += (p - phase_resolution)**2

        fitness = -100 * count / len(freq)
        return (fitness,)

    except Exception as e:
        print(f"Error in run_program: {e}")
        return (-1e10,)

# ============================================================
# DEAP GA INITIALIZATION AND PARAMETERS
# ============================================================

seq1 = [4, 5]

if not hasattr(creator, "FitnessMax"):
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
if not hasattr(creator, "Individual"):
    creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("attr_val", random.choice, seq1)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_val, 40)
toolbox.register("population", tools.initRepeat, list, toolbox.individual, MU)
toolbox.register("evaluate", run_program)
toolbox.register("mate", tools.cxUniform, indpb=0.15)
toolbox.register("mutate", tools.mutUniformInt, low=4, up=5, indpb=0.15)
toolbox.register("select", tools.selBest)

stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", np.mean)
stats.register("min", np.min)
stats.register("max", np.max)

LOGFILE = "ga_log.csv"

# ============================================================
# MAIN GA LOOP
# ============================================================

def main():

    #start with gen=0 with option to seed or not then continue the rest
    start_gen = 0
    population = None
    hof = tools.HallOfFame(1)
    with open(LOGFILE, 'w') as f:
        f.write(f"{'Gen':>4} | {'Avg Fit':>10} | {'Min Fit':>10} | {'Max Fit':>10} | {'Best Overall':>13}\n")
        f.write("-" * 60 + "\n")

    if SEED:
        print("Evaluating initial population with seeding...")
        seed_file = "seeds.npy"
        if os.path.exists(seed_file):
            print(f"Loading population from {seed_file}")
            pop_data = np.load(seed_file, allow_pickle=True)
            population = [creator.Individual([int(x) for x in ind]) for ind in pop_data]

    if population is None:
        population = toolbox.population()
        print("Evaluating initial population...")

    gen_folder = "generation_init"
    os.makedirs(gen_folder, exist_ok=True)

    np.save(os.path.join(gen_folder, "offspring.npy"), population, allow_pickle=True)

    cmd = ["sbatch", f"--array=0-{MU-1}", "run_generation.sh", gen_folder]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    print("Submitted initial population:", result.stdout.strip())

    start_wait = time.time()
    while True:
        done = all(os.path.exists(os.path.join(gen_folder, f"fit_{i}.npy")) for i in range(MU))
        if done:
            break
        if time.time() - start_wait > 6 * 3600:
            raise TimeoutError("Timed out waiting for initial population results.")
        time.sleep(CHECK_INTERVAL)

    for i, ind in enumerate(population):
        fit = np.load(os.path.join(gen_folder, f"fit_{i}.npy"))[0]
        ind.fitness.values = (fit,)

    if AUTO_CLEANUP:
        shutil.rmtree(gen_folder, ignore_errors=True)

    start_gen = 1
    hof.update(population)

    #now the rest of them after initialization
    for gen in range(start_gen, NGEN):
        print(f"\n{'='*60}")
        print(f"Generation {gen}/{NGEN-1}")
        print(f"{'='*60}")

        gen_folder = f"generation_{gen}"
        os.makedirs(gen_folder, exist_ok=True)

        offspring = algorithms.varAnd(population, toolbox, CXPB, MUTPB)
        np.save(os.path.join(gen_folder, "offspring.npy"), offspring, allow_pickle=True)
        print(f"Generated {len(offspring)} offspring")

        cmd = ["sbatch", f"--array=0-{len(offspring)-1}", "run_generation.sh", gen_folder]
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            job_id = result.stdout.strip().split()[-1]
            print(f"Submitted job: {job_id}")
        except subprocess.CalledProcessError as e:
            print(f"SLURM submission failed: {e.stderr}")
            print("Skipping this generation to avoid resubmission loop.")
            continue

        print("Waiting for fitness results...")
        start_wait = time.time()
        last_count = 0

        while True:
            completed = sum(1 for i in range(len(offspring))
                            if os.path.exists(os.path.join(gen_folder, f"fit_{i}.npy")))
            if completed > last_count:
                print(f"  Progress: {completed}/{len(offspring)} completed")
                last_count = completed

            if completed == len(offspring):
                print("All fitness evaluations complete!")
                break

            elapsed = time.time() - start_wait
            if elapsed > 6 * 3600:
                raise TimeoutError(f"Timed out waiting for {gen_folder} results after {elapsed/3600:.1f} hours.")
            time.sleep(CHECK_INTERVAL)

        fits = []
        for i in range(len(offspring)):
            fit_file = os.path.join(gen_folder, f"fit_{i}.npy")
            fit = np.load(fit_file)[0]
            fits.append(fit)

        for ind, fit in zip(offspring, fits):
            ind.fitness.values = (fit,)

        population[:] = toolbox.select(population + offspring, MU)
        hof.update(population)

        # Remove old population and hof
        for old_ckpt in [f for f in os.listdir(".") if (f.startswith("population_gen_") or f.startswith("hof_gen_")) and f.endswith(".npy")]:
            os.remove(old_ckpt)

        # Save new population and hof
        np.save(f"population_gen_{gen}.npy", [list(ind) for ind in population], allow_pickle=True)
        np.save(f"hof_gen_{gen}.npy", [list(ind) for ind in hof], allow_pickle=True)

        record = stats.compile(population)
        with open(LOGFILE, 'a') as f:
            f.write(f"{gen:4d} | {record['avg']:10.3f} | {record['min']:10.3f} | {record['max']:10.3f} | {hof[0].fitness.values[0]:13.3f}\n")


        print(f"\nGeneration {gen} Statistics:")
        print(f"  Average fitness: {record['avg']:.3f}")
        print(f"  Min fitness:     {record['min']:.3f}")
        print(f"  Max fitness:     {record['max']:.3f}")
        print(f"  Best overall:    {hof[0].fitness.values[0]:.3f}")
        #print(f"  Best params:     {list(hof[0])}")

        if AUTO_CLEANUP:
            print(f"Cleaning up {gen_folder}...")
            shutil.rmtree(gen_folder, ignore_errors=True)
            print(f"Cleaning up ga_eval folder of slurm files...")
            shutil.rmtree('ga_eval', ignore_errors=True)

    print("\n" + "="*60)
    print("OPTIMIZATION COMPLETE")
    print("="*60)
    print(f"Best fitness: {hof[0].fitness.values[0]:.3f}")
    print(f"Best individual: {list(hof[0])}")

if __name__ == "__main__":
    main()
