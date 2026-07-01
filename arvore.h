#ifndef ARVORE_H
#define ARVORE_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "geometria.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Estrutura do Nó expandida com propriedades físicas (Seção 3)
typedef struct No {
    struct No *esq;
    struct No *dir;
    struct No *pai;
    int id;
    Point p;               // Coordenada distal do segmento
    
    // Propriedades físicas do segmento (pai -> p)
    double raio;           // Raio do segmento
    double comprimento;    // Comprimento do segmento
    double fluxo;          // Fluxo no segmento
    double resistencia;    // Resistência hidráulica (Poiseuille)
    double volume;         // Volume intravascular do segmento
    int qtd_term_distal;   // Quantidade de terminais distais
} No;

typedef No* ptrNo;

// Funções base de gerenciamento da estrutura (MiniCCO-0)
ptrNo  cria_no(Point p, ptrNo pai, int id);
bool   insere_filho(ptrNo pai, ptrNo filho);
void   insere_bifurcacao(ptrNo pai, ptrNo filho_antigo, ptrNo bifurc, ptrNo novo_terminal);
void   percurso_nos(ptrNo raiz, void (*callback)(ptrNo));
int    conta_folhas(ptrNo raiz);
double comprimento_total(ptrNo raiz);
void   libera_arvore(ptrNo raiz);

// Novas funções físicas e hidráulicas (MiniCCO-1)
double calculaComprimento(ptrNo seg);
double calculaResistencia(double mu, double comprimento, double raio);
double calculaVolume(double comprimento, double raio);
double calculaVolumeTotal(ptrNo raiz);

int    atualizaQtdTerminaisDistais(ptrNo no);
void   atualizaFluxos(ptrNo no, double Qterm);
void   atualizaRaiosPorFluxo(ptrNo no, double gamma);
void   atualizaGeometriaFisica(ptrNo raiz, double Qterm, double gamma, double mu);
double funcaoCustoVolume(ptrNo raiz);

#endif
