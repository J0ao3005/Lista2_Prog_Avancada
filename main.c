#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "arvore.h"
#include "candidatos.h"
#include "exportar.h" 

// Função auxiliar para calcular o raio médio dos segmentos
static void soma_raios_segmentos(ptrNo no, double *soma, int *qtd) {
    if (no == NULL) return;
    if (no->pai != NULL) {
        *soma += no->raio;
        (*qtd)++;
    }
    soma_raios_segmentos(no->esq, soma, qtd);
    soma_raios_segmentos(no->dir, soma, qtd);
}

int main(int argc, char *argv[]) {

    // ── Leitura dos novos argumentos (Seção 13) ─────────────────────
    if (argc != 5) {
        fprintf(stderr, "Uso: %s <Nterm> <R> <gamma> <M>\n", argv[0]);
        return 1;
    }

    int    Nterm   = atoi(argv[1]);
    double R       = atof(argv[2]);
    double gamma   = atof(argv[3]);
    int    M       = atoi(argv[4]);
    double epsilon = 0.01 * R;

    // Constantes físicas (Seção 4)
    double Q_perf = 8.33e-6; 
    double mu     = 3.6e-3;
    double Q_term = Q_perf / Nterm;

    printf("=== MiniCCO-1 ===\n");
    printf("Nterm   : %d\n", Nterm);
    printf("R       : %.4f\n", R);
    printf("Gamma   : %.2f\n", gamma);
    printf("M (res) : %d\n\n", M);

    // ── Inicialização ───────────────────────────────────────────────
    srand((unsigned int)time(NULL));
    clock_t inicio_tempo = clock();

    Point raiz_pt  = {0.0, R}; 
    ptrNo raiz     = cria_no(raiz_pt, NULL, 0);
    int   proximo_id = 1;
    int   rejeitados = 0;
    int   inseridos  = 0;

    // ── Loop principal de crescimento ───────────────────────────────
    for (int i = 0; i < Nterm; i++) {
        Point novo = gera_ponto_aleatorio(R);

        // Caso especial: árvore só tem a raiz
        if (raiz->esq == NULL && raiz->dir == NULL) {
            ptrNo novo_no = cria_no(novo, raiz, proximo_id++);
            if (novo_no == NULL) continue;

            insere_filho(raiz, novo_no);
            
            // Importante: Atualizar física já no primeiro segmento
            atualizaGeometriaFisica(raiz, Q_term, gamma, mu);
            inseridos++;
            continue;
        }

        // Caso normal: encontra o melhor segmento com otimização baricêntrica
        Point melhor_X;
        ptrNo filho = seleciona_melhor_cco1(novo, raiz, R, epsilon, M, 
                                            Q_term, gamma, mu, 
                                            &rejeitados, &melhor_X);

        if (filho == NULL) continue;

        // Efetiva a inserção no ponto ótimo (melhor_X) encontrado na grade
        ptrNo bifurc  = cria_no(melhor_X, filho->pai, proximo_id++);
        ptrNo novo_no = cria_no(novo, bifurc, proximo_id++);

        insere_bifurcacao(filho->pai, filho, bifurc, novo_no);
        
        // Atualiza a física e os volumes da árvore para o próximo ciclo
        atualizaGeometriaFisica(raiz, Q_term, gamma, mu);
        inseridos++;
    }

    clock_t fim_tempo = clock();
    double tempo_execucao = (double)(fim_tempo - inicio_tempo) / CLOCKS_PER_SEC;

    // ── Estatísticas finais obrigatórias (Seção 14) ─────────────────
    double soma_raios = 0.0;
    int qtd_segmentos = 0;
    soma_raios_segmentos(raiz, &soma_raios, &qtd_segmentos);
    double raio_medio = (qtd_segmentos > 0) ? (soma_raios / qtd_segmentos) : 0.0;

    // O raio da raiz é representado pelo raio do segmento que sai do nó raiz
    double raio_raiz = 0.0;
    if (raiz->esq) raio_raiz = raiz->esq->raio;

    printf("=== Resultados ===\n");
    printf("Nos totais          : %d\n", proximo_id);
    printf("Segmentos totais    : %d\n", qtd_segmentos);
    printf("Terminais inseridos : %d\n", inseridos);
    printf("Comprimento total   : %.6f\n", comprimento_total(raiz));
    printf("Volume intravascular: %.10e\n", funcaoCustoVolume(raiz));
    printf("Raio da raiz        : %.10e\n", raio_raiz);
    printf("Raio medio          : %.10e\n", raio_medio);
    printf("Conexoes rejeitadas : %d\n", rejeitados);
    printf("Tempo de execucao   : %.3f segundos\n\n", tempo_execucao);

    // ── Exportação ──────────────────────────────────────────────────
    exporta_csv_fisica(raiz, "arvore_fisica.csv");

    // ── Liberação de memória ────────────────────────────────────────
    libera_arvore(raiz);

    return 0;
}
