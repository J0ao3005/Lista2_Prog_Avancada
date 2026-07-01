import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else "arvore_fisica.csv"
    
    try:
        dados = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Erro: arquivo '{filename}' nao encontrado.")
        sys.exit(1)

    if dados.empty:
        print("Nenhum segmento encontrado no arquivo.")
        sys.exit(1)

    print(f"Lendo {len(dados)} segmentos de {filename}...")

    plt.figure(figsize=(10, 10))
    
    # Encontra o maior e o menor raio para normalizar a espessura da linha
    raio_max = dados['raio'].max()
    raio_min = dados['raio'].min()
    
    # Plota segmento por segmento
    for index, row in dados.iterrows():
        x = [row['x0'], row['x1']]
        y = [row['y0'], row['y1']]
        raio = row['raio']
        
        # Normalizando a espessura da linha entre 1.0 (mais fina) e 8.0 (mais grossa)
        if raio_max > 0 and raio_max != raio_min:
            espessura = 1.0 + 7.0 * ((raio - raio_min) / (raio_max - raio_min))
        else:
            espessura = 3.0 # Fator padrão caso a árvore só tenha um segmento
            
        # Desenha a "artéria" (linha com espessura variável)
        plt.plot(x, y, color='crimson', linewidth=espessura, zorder=1, solid_capstyle='round')
        
        # Desenha o nó na ponta
        plt.scatter([row['x1']], [row['y1']], color='steelblue', s=15, zorder=2)

    # Destaca a raiz (usa a coordenada x0, y0 do primeiro segmento)
    raiz_x = dados.iloc[0]['x0']
    raiz_y = dados.iloc[0]['y0']
    plt.scatter([raiz_x], [raiz_y], color='gold', edgecolor='black', s=150, zorder=3, label='Raiz')

    plt.title("MiniCCO-1 — Árvore com Escala de Raios")
    plt.axis('equal')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    
    # Salva em arquivo PNG em vez de tentar abrir a interface gráfica
    plt.savefig("arvore_fisica_vis.png", dpi=300, bbox_inches='tight')
    print("Sucesso! Imagem salva como 'arvore_fisica_vis.png'.")

if __name__ == "__main__":
    main()
