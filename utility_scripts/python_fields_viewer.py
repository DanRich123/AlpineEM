import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap
import time as tp
import os
import subprocess

# SEE LICENSE FILE FOR LICENSE INFORMATION

# This script generates videos or images of the fields saved by the user for the 2D slice field data output option
# User only needs to modify the setup section.
# E,H fields are cell centered and at the same time.

###############################################################################
#### SETUP ####################################################################
Type_sim = 'time domain'        # 'time domain' or 'frequency domain'
target_freq = 5E9               # target frequency if frequency domain, else unused
still_frame = True              # still frame if frequency domain, else unused
data_filename = 'data.dat'      # file to pull info from

Fields = ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz']
scale = 'auto'                  # 'auto' or 'manual'
# If manual, these bounds are used:
max_vals = [1, 1, 1, 1/377, 1/377, 1/377]
min_vals = [-1, -1, -1, -1/377, -1/377, -1/377]

color_mapping = 'nipy_spectral' # color mapping for the fields
Type_fields = 'total'           # 'total', 'incident', 'scattered'
clear_filename = 'clear.dat'    # clear geometry file name to pule info from

use_geometry = True    # show geometry of cross section or not
opacity = 0.4           # opacity for the combination: 0.0 shows no EM fields, 1.0 shows no geometry
geom_colors = "Pastel1" # color mapping for the geometry
#colors = ["gray", "black"]
#geom_colors = LinearSegmentedColormap.from_list("custom_div", colors)

KMAX = True
###############################################################################
###############################################################################

###############################################################################
# Determine output writer: prefer ffmpeg (MP4), fall back to pillow (GIF)
def get_writer(still):
    if still:
        return None, '.png'
    if animation.FFMpegWriter.isAvailable():
        return animation.FFMpegWriter(fps=30), '.mp4'
    return animation.PillowWriter(fps=30), '.gif'

#setup ffmpeg settings
def make_ffmpeg_proc(out_path, fig, fps=30):
    w = int(fig.get_figwidth()  * fig.dpi)
    h = int(fig.get_figheight() * fig.dpi)
    # round to even (ffmpeg requirement)
    w += w % 2
    h += h % 2
    cmd = [
    'ffmpeg', '-y',
    '-f', 'rawvideo',
    '-vcodec', 'rawvideo',
    '-s', f'{w}x{h}',
    '-pix_fmt', 'rgba',
    '-r', str(fps),
    '-i', 'pipe:0',
    '-vcodec', 'libopenh264',
    '-pix_fmt', 'yuv420p',
    out_path
    ]
    return subprocess.Popen(cmd, stdin=subprocess.PIPE), w, h

###############################################################################
# Load header constants
start_time = tp.time()

x_size        = int(np.loadtxt(data_filename, usecols=3, skiprows=0, max_rows=1))
y_size        = int(np.loadtxt(data_filename, usecols=3, skiprows=1, max_rows=1))
z_size        = int(np.loadtxt(data_filename, usecols=3, skiprows=2, max_rows=1))
dx            = np.loadtxt(data_filename, usecols=3, skiprows=3, max_rows=1)
dy            = np.loadtxt(data_filename, usecols=4, skiprows=3, max_rows=1)
dz            = np.loadtxt(data_filename, usecols=5, skiprows=3, max_rows=1)
del_t         = np.loadtxt(data_filename, usecols=3, skiprows=4, max_rows=1)
time_blocks   = int(np.loadtxt(data_filename, usecols=5, skiprows=5, max_rows=1))
slice_type    = int(np.loadtxt(data_filename, usecols=2, skiprows=9, max_rows=1))
slice_location= np.loadtxt(data_filename, usecols=3, skiprows=10, max_rows=1)
video1        = int(np.loadtxt(data_filename, usecols=3, skiprows=11, max_rows=1))
video2        = int(np.loadtxt(data_filename, usecols=3, skiprows=12, max_rows=1))

###############################################################################
# Axis labels and limits based on slice type
slice_meta = {
    0: ('Cubes in Y direction', 'Cubes in Z direction', 'X',
        (0.5, y_size + 0.5), (0.5, z_size + 0.5)),
    1: ('Cubes in X direction', 'Cubes in Z direction', 'Y',
        (0.5, x_size + 0.5), (0.5, z_size + 0.5)),
    2: ('Cubes in X direction', 'Cubes in Y direction', 'Z',
        (0.5, x_size + 0.5), (0.5, y_size + 0.5)),
}
xlabel, ylabel, zlabel, xlim, ylim = slice_meta[slice_type]

###############################################################################
# Load geometry once (shared across all field components)
geom_rgba = None
if use_geometry:
    geom = np.fromfile('{}'.format(data_filename[:-4] + '_geometry.bin'), dtype=np.float32)
    geom = geom.reshape((4, x_size, y_size, z_size), order='F')
    materials_grid, sheetx, sheety, sheetz = [np.ascontiguousarray(a) for a in geom]

    sl = int(slice_location) - 1
    if slice_type == 0:
        mat    = materials_grid[sl, :, :]
        sheet1 = sheety[sl, :, :]
        sheet2 = sheetz[sl, :, :]
        dims   = (y_size, z_size)
        s1_dir, s2_dir = 'y', 'z'
    elif slice_type == 1:
        mat    = materials_grid[:, sl, :]
        sheet1 = sheetx[:, sl, :]
        sheet2 = sheetz[:, sl, :]
        dims   = (x_size, z_size)
        s1_dir, s2_dir = 'x', 'z'
    else:
        mat    = materials_grid[:, :, sl]
        sheet1 = sheetx[:, :, sl]
        sheet2 = sheety[:, :, sl]
        dims   = (x_size, y_size)
        s1_dir, s2_dir = 'x', 'y'

    # Build geometry lines list once
    geom_lines = []
    for arr, axis in [(sheet1, 'first'), (sheet2, 'second')]:
        active_rows = np.where(np.any(arr != 0, axis=1))[0]
        for ii in active_rows:
            active_cols = np.where(arr[ii, :] != 0)[0]
            for jj in active_cols:
                if axis == 'first':
                    geom_lines.append(([ii+0.5, ii+0.5], [jj+0.5, jj+1.5]))
                else:
                    geom_lines.append(([ii+0.5, ii+1.5], [jj+0.5, jj+0.5]))

###############################################################################
# Auto-scale pass: load all data first, compute bounds, then animate
# This avoids a second file read when scale='auto'.

if scale == 'auto':
    min_vals = []
    max_vals = []

    for field in Fields:
        file_path = '{}_video_{}.bin'.format(data_filename[:-4], field)
        if KMAX:
            d = np.fromfile(file_path, dtype=np.complex64)[1:]
        else:
            d = np.fromfile(file_path, dtype=np.float32)[1:-1]
        d = d.reshape(time_blocks, video2, video1)

        if Type_fields in ('scatter', 'incident'):
            cp = '{}_video_{}.bin'.format(clear_filename[:-4], field)
            if KMAX:
                c = np.fromfile(cp, dtype=np.complex64)[1:]
            else:
                c = np.fromfile(cp, dtype=np.float32)[1:-1]
            c = c.reshape(time_blocks, video2, video1)
            d = c if Type_fields == 'incident' else d - c

        if Type_sim == 'time domain':
            min_vals.append(float(np.real(d).min()))
            max_vals.append(float(np.real(d).max()))
        else:
            dfft = np.fft.fft(d, axis=0)
            freq = np.fft.fftfreq(time_blocks, del_t)
            idx  = np.argmin(np.abs(freq - target_freq))
            min_vals.append(float(np.real(dfft[idx]).min()))
            max_vals.append(float(np.real(dfft[idx]).max()))

###############################################################################
# Main loop: one figure per field component
writer, ext = get_writer(still_frame if Type_sim == 'frequency domain' else False)
# For time domain we never want a still PNG
if Type_sim == 'time domain':
    writer, ext = get_writer(False)

for i, field in enumerate(Fields):

    # --- Load data ---
    file_path = '{}_video_{}.bin'.format(data_filename[:-4], field)
    if KMAX:
        data = np.fromfile(file_path, dtype=np.complex64)[1:]
    else:
        data = np.fromfile(file_path, dtype=np.float32)[1:-1]
    data = data.reshape(time_blocks, video2, video1)

    if Type_fields in ('scatter', 'incident'):
        cp = '{}_video_{}.bin'.format(clear_filename[:-4], field)
        if KMAX:
            clear = np.fromfile(cp, dtype=np.complex64)[1:]
        else:
            clear = np.fromfile(cp, dtype=np.float32)[1:-1]
        clear = clear.reshape(time_blocks, video2, video1)
        data = clear if Type_fields == 'incident' else data - clear

    vmin, vmax = min_vals[i], max_vals[i]

    # --- Pre-compute all frame arrays as one stacked array ---
    if Type_sim == 'time domain':
        # shape: (time_blocks, video2, video1)
        all_frames = np.real(data)
        n_frames   = time_blocks

    else:  # frequency domain
        dfft  = np.fft.fft(data, axis=0)
        freq  = np.fft.fftfreq(time_blocks, del_t)
        idx   = np.argmin(np.abs(freq - target_freq))
        n_frames = 1 if still_frame else 100
        phases   = np.linspace(0, 2 * np.pi, n_frames)
        # broadcast: (n_frames, video2, video1)
        all_frames = np.real(dfft[idx][np.newaxis, :, :] *
                             np.exp(1j * phases)[:, np.newaxis, np.newaxis])

    # --- Build figure ---
    fig, ax = plt.subplots()

    if use_geometry:
        geom_im = ax.imshow(mat.T,
                            origin='lower',
                            extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
                            aspect='equal',
                            cmap=geom_colors,
                            interpolation='nearest',
                            zorder=0)

        for (xs, ys) in geom_lines:
            ax.plot(xs, ys, color='black', linewidth=0.8, zorder=2)

    im = ax.imshow(all_frames[0],
                origin='lower',
                extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
                aspect='equal',
                cmap=color_mapping,
                vmin=vmin, vmax=vmax,
                interpolation='nearest',
                alpha=opacity if use_geometry else 1.0,
                zorder=1)

    cb_label = 'V/m' if field in ('Ex', 'Ey', 'Ez') else 'A/m'
    fig.colorbar(im, ax=ax, label=cb_label)
    fig.set_layout_engine('none') 

    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.grid(True, color='lightgray', which='both', linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)

    if KMAX and Type_sim == 'time domain':
        title = 'Real amplitude {} - cube {} - {} direction'.format(
            field, int(slice_location), zlabel)
    else:
        title = 'Amplitude {} - cube {} - {} direction'.format(
            field, int(slice_location), zlabel)

    ax.set_title(title, x=0.5, y=1.05, ha='center')

    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(ax.bbox)

    if ext == '.png':
        fig.savefig('{}{}'.format(field, ext), dpi=150)
        plt.close(fig)
        continue

    proc, fw, fh = make_ffmpeg_proc('{}{}'.format(field, ext), fig, fps=30)

    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(ax.bbox)

    for frame_idx in range(n_frames):
        fig.canvas.restore_region(background)
        im.set_data(all_frames[frame_idx])
        ax.draw_artist(im)
        fig.canvas.blit(ax.bbox)
        buf = fig.canvas.buffer_rgba()
        proc.stdin.write(buf.tobytes())

    proc.stdin.close()
    proc.wait()
    plt.close(fig)

    if i < len(Fields) - 1:
        print('{} done. Starting {}.'.format(field, Fields[i + 1]))

###############################################################################
end_time = tp.time()
print('Total time: {:.2f} minutes'.format((end_time - start_time) / 60))