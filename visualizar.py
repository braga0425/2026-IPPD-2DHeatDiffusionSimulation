import glob
import os
import re
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")


COLORMAP = "inferno"


def ler_instante(caminho):
    with open(caminho) as f:
        cab = f.readline()
    m = re.search(r"passo=(\d+)\s+tempo=([\d.eE+-]+)", cab)
    passo = int(m.group(1)) if m else -1
    tempo = float(m.group(2)) if m else float("nan")
    matriz = np.loadtxt(caminho)
    return passo, tempo, matriz


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dir_saida = args[0] if args else "saida"
    fazer_gif = "--gif" in sys.argv

    arquivos = sorted(glob.glob(os.path.join(dir_saida, "passo_*.txt")))
    if not arquivos:
        sys.exit(f"Nenhum arquivo passo_*.txt encontrado em '{dir_saida}'.")

    instantes = [ler_instante(a) for a in arquivos]

    vmin = min(m.min() for _, _, m in instantes)
    vmax = max(m.max() for _, _, m in instantes)

    dir_img = os.path.join(dir_saida, "imagens")
    os.makedirs(dir_img, exist_ok=True)

    quadros = []
    for passo, tempo, matriz in instantes:
        fig, ax = plt.subplots(figsize=(6.4, 5.4))
        im = ax.imshow(matriz, cmap=COLORMAP, vmin=vmin, vmax=vmax,
                       origin="upper", interpolation="nearest")
        ax.set_title(f"Difusão de calor 2D — passo {passo} (t = {tempo:.4f} s)")
        ax.set_xlabel("j (coluna)")
        ax.set_ylabel("i (linha)")
        fig.colorbar(im, ax=ax, label="Temperatura")
        fig.tight_layout()

        png = os.path.join(dir_img, f"passo_{passo:06d}.png")
        fig.savefig(png, dpi=120)
        plt.close(fig)
        quadros.append((passo, tempo, matriz))
        print(f"gerado {png}")

    if fazer_gif:
        from matplotlib.animation import FuncAnimation, PillowWriter

        fig, ax = plt.subplots(figsize=(6.4, 5.4))
        im = ax.imshow(quadros[0][2], cmap=COLORMAP, vmin=vmin, vmax=vmax,
                       origin="upper", interpolation="nearest")
        fig.colorbar(im, ax=ax, label="Temperatura")
        ax.set_xlabel("j (coluna)")
        ax.set_ylabel("i (linha)")

        def atualizar(k):
            passo, tempo, matriz = quadros[k]
            im.set_data(matriz)
            ax.set_title(f"Difusão de calor 2D — passo {passo} (t = {tempo:.4f} s)")
            return (im,)

        anim = FuncAnimation(fig, atualizar, frames=len(quadros))
        gif = os.path.join(dir_saida, "animacao.gif")
        anim.save(gif, writer=PillowWriter(fps=5))
        plt.close(fig)
        print(f"gerado {gif}")


if __name__ == "__main__":
    main()
