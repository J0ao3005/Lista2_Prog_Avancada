#include "arvore.h"

ptrNo cria_no(Point p, ptrNo pai, int id) {
    No *a = malloc(sizeof(No));
    if (a == NULL) {
        fprintf(stderr, "Erro: falha ao alocar no.\n");
        return NULL;
    }
    a->esq = NULL;
    a->dir = NULL;
    a->pai = pai;
    a->p   = p;
    a->id  = id;
    
    // Inicialização dos campos físicos para o MiniCCO-1
    a->raio = 0.0;
    a->comprimento = 0.0;
    a->fluxo = 0.0;
    a->resistencia = 0.0;
    a->volume = 0.0;
    a->qtd_term_distal = 0;
    
    return a;
}

bool insere_filho(ptrNo pai, ptrNo filho) {
    if (pai->esq == NULL) { pai->esq = filho; return true; }
    if (pai->dir == NULL) { pai->dir = filho; return true; }
    return false;
}

void insere_bifurcacao(ptrNo pai, ptrNo filho_antigo, ptrNo bifurc, ptrNo novo_terminal) {
    if (pai->esq == filho_antigo)
        pai->esq = bifurc;
    else
        pai->dir = bifurc;

    bifurc->pai = pai;
    bifurc->esq = filho_antigo;
    bifurc->dir = novo_terminal;

    filho_antigo->pai  = bifurc;
    novo_terminal->pai = bifurc;
}

void percurso_nos(ptrNo raiz, void (*callback)(ptrNo)) {
    if (raiz == NULL) return;
    callback(raiz);
    percurso_nos(raiz->esq, callback);
    percurso_nos(raiz->dir, callback);
}

int conta_folhas(ptrNo raiz) {
    if (raiz == NULL) return 0;
    if (raiz->esq == NULL && raiz->dir == NULL) return 1;
    return conta_folhas(raiz->esq) + conta_folhas(raiz->dir);
}

double comprimento_total(ptrNo raiz) {
    if (raiz == NULL) return 0;
    double comp = 0;
    if (raiz->pai != NULL)
        comp = distancia(raiz->p, raiz->pai->p);
    return comp + comprimento_total(raiz->esq) + comprimento_total(raiz->dir);
}

void libera_arvore(ptrNo raiz) {
    if (raiz == NULL) return;
    libera_arvore(raiz->esq);
    libera_arvore(raiz->dir);
    free(raiz);
}

// ─────────────────────────────────────────────────────────────────────
// Funções Físicas e Hidráulicas (MiniCCO-1)
// ─────────────────────────────────────────────────────────────────────

// Parte A: Cálculos geométricos e físicos estruturais
double calculaComprimento(ptrNo seg) {
    if (seg == NULL || seg->pai == NULL) return 0.0;
    return distancia(seg->p, seg->pai->p);
}

double calculaResistencia(double mu, double comprimento, double raio) {
    if (raio <= 0.0) return 0.0;
    // Lei de Poiseuille: R = (8 * mu * L) / (pi * r^4)
    return (8.0 * mu * comprimento) / (M_PI * pow(raio, 4.0));
}

double calculaVolume(double comprimento, double raio) {
    // V = pi * r^2 * L
    return M_PI * raio * raio * comprimento;
}

double calculaVolumeTotal(ptrNo raiz) {
    if (raiz == NULL) return 0.0;
    double vol = 0.0;
    if (raiz->pai != NULL) {
        vol = raiz->volume;
    }
    return vol + calculaVolumeTotal(raiz->esq) + calculaVolumeTotal(raiz->dir);
}

// Parte B: Contagem de Terminais Distais (Percurso pós-ordem)
int atualizaQtdTerminaisDistais(ptrNo no) {
    if (no == NULL) return 0;
    
    if (no->esq == NULL && no->dir == NULL) {
        no->qtd_term_distal = 1;
    } else {
        no->qtd_term_distal = atualizaQtdTerminaisDistais(no->esq) + 
                              atualizaQtdTerminaisDistais(no->dir);
    }
    return no->qtd_term_distal;
}

// Parte C: Fluxo em Cada Segmento
void updateFluxosRecursivo(ptrNo no, double Qterm) {
    if (no == NULL) return;
    no->fluxo = no->qtd_term_distal * Qterm;
    updateFluxosRecursivo(no->esq, Qterm);
    updateFluxosRecursivo(no->dir, Qterm);
}

void atualizaFluxos(ptrNo no, double Qterm) {
    updateFluxosRecursivo(no, Qterm);
}

// Parte D: Lei de Bifurcação (Escala dos raios com C = 1)
void updateRaiosRecursivo(ptrNo no, double gamma) {
    if (no == NULL) return;
    if (gamma > 0.0 && no->fluxo > 0.0) {
        no->raio = pow(no->fluxo, 1.0 / gamma);
    } else {
        no->raio = 0.0;
    }
    updateRaiosRecursivo(no->esq, gamma);
    updateRaiosRecursivo(no->dir, gamma);
}

void atualizaRaiosPorFluxo(ptrNo no, double gamma) {
    updateRaiosRecursivo(no, gamma);
}

// Função auxiliar para recalcular as demais propriedades dependentes do raio
static void atualizaPropriedadesSegmento(ptrNo no, double mu) {
    if (no == NULL) return;
    
    no->comprimento = calculaComprimento(no);
    no->resistencia = calculaResistencia(mu, no->comprimento, no->raio);
    no->volume      = calculaVolume(no->comprimento, no->raio);
    
    atualizaPropriedadesSegmento(no->esq, mu);
    atualizaPropriedadesSegmento(no->dir, mu);
}

// Consolidação do pipeline de atualização geométrica e física
void atualizaGeometriaFisica(ptrNo raiz, double Qterm, double gamma, double mu) {
    if (raiz == NULL) return;
    
    atualizaQtdTerminaisDistais(raiz);
    atualizaFluxos(raiz, Qterm);
    atualizaRaiosPorFluxo(raiz, gamma);
    atualizaPropriedadesSegmento(raiz, mu);
}

// Parte E: Função Custo baseada em Volume Intravascular
double funcaoCustoVolume(ptrNo raiz) {
    return calculaVolumeTotal(raiz);
}
