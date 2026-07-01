import sys
import numpy as np
import pandas as pd
import pyvista as pv

def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else "arvore_fisica.csv"
    
    try:
        dados = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Erro: arquivo '{filename}' nao encontrado.")
        print("Execute primeiro: ./minicco1 <Nterm> <R> <gamma> <M>")
        sys.exit(1)
        
    if dados.empty:
        print("Nenhum segmento encontrado no arquivo.")
        sys.exit(1)

    print(f"Lendo {len(dados)} segmentos de {filename}...")

    # Configuração do PyVista
    plotter = pv.Plotter()
    plotter.background_color = 'white'
    
    pontos_nos = []
    raiz_pt = None

    # Monta os tubos cilíndricos segmento por segmento
    for index, row in dados.iterrows():
        p0 = np.array([row["x0"], row["y0"], 0.0])
        p1 = np.array([row["x1"], row["y1"], 0.0])
        raio = row["raio"]
        
        # Garante que raios extremamente pequenos não "sumam" na visualização
        if raio < 1e-6:
            raio = 1e-6
            
        # Cria a linha e aplica a espessura (tubo)
        linha = pv.Line(p0, p1)
        tubo = linha.tube(radius=raio, n_sides=20)
        
        plotter.add_mesh(tubo, color='crimson', smooth_shading=True)
        
        # Armazena os pontos para desenhar as esferas depois
        pontos_nos.append(p1)
        if raiz_pt is None:
            raiz_pt = p0  # O primeiro x0, y0 do arquivo será a raiz
            pontos_nos.append(p0)

    # Adiciona as esferas nos nós e bifurcações
    if pontos_nos:
        nuvem_pontos = pv.PolyData(np.array(pontos_nos))
        plotter.add_mesh(nuvem_pontos, point_size=5, render_points_as_spheres=True, color='steelblue', label='Nós/Bifurcações')

    # Destaca a raiz
    if raiz_pt is not None:
        plotter.add_mesh(pv.PolyData(np.array([raiz_pt])), point_size=15, render_points_as_spheres=True, color='gold', label='Raiz')

    # Configura a câmera e exibe
    plotter.view_xy()
    plotter.add_legend()
    plotter.show(title="MiniCCO-1 — Árvore com Escala de Raios")

if __name__ == "__main__":
    main()
