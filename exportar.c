#include "exportar.h"
#include <stdio.h>

static FILE *arq_saida;

static void escreve_segmento_fisico(ptrNo no) {
    if (no == NULL || no->pai == NULL) return; // A raiz em si não forma segmento

    int id_pai = no->pai->id;
    int id_no  = no->id;

    fprintf(arq_saida, "%d,%d,%.6f,%.6f,%.6f,%.6f,%.10e,%.6f,%.10e,%.10e,%.10e\n",
            id_no, id_pai,
            no->pai->p.x, no->pai->p.y,
            no->p.x, no->p.y,
            no->raio,
            no->comprimento,
            no->fluxo,
            no->resistencia,
            no->volume);
}

void exporta_csv_fisica(ptrNo raiz, const char *filename) {
    arq_saida = fopen(filename, "w");
    if (arq_saida == NULL) {
        fprintf(stderr, "Erro: nao foi possivel abrir '%s' para escrita.\n", filename);
        return;
    }

    // Cabeçalho exigido pelo MiniCCO-1
    fprintf(arq_saida, "id,pai,x0,y0,x1,y1,raio,comprimento,fluxo,resistencia,volume\n");

    percurso_nos(raiz, escreve_segmento_fisico);

    fclose(arq_saida);
    printf("Arvore exportada com fisica para '%s'\n", filename);
}
