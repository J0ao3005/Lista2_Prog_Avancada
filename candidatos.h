#ifndef CANDIDATOS_H
#define CANDIDATOS_H

#include "arvore.h"
#include <float.h>

#define MAX_NOS 10000

// Gera ponto aleatório dentro do domínio circular
Point gera_ponto_aleatorio(double R);

// Coleta todos os nós não-raiz da árvore
void coleta_filhos(ptrNo raiz, ptrNo *lista, int *total);

// Coleta todos os segmentos existentes para checagem de colisão
void coleta_segmentos(ptrNo raiz, Segment *lista, int *total);

// ─────────────────────────────────────────────────────────────────────
// Funções Matemáticas do MiniCCO-1 (Busca em Grade Baricêntrica)
// ─────────────────────────────────────────────────────────────────────

// Calcula um ponto X usando coordenadas baricêntricas
Point pontoBaricentrico(Point A, Point B, Point C, double alpha, double beta, double lambda);

// Verifica se P está dentro do triângulo ABC usando produto vetorial (orientação)
int pontoDentroTriangulo(Point A, Point B, Point C, Point P);

// Percorre a grade de resolução M no triângulo, avaliando o volume total
Point otimizaBifurcacaoPorGrade(ptrNo raiz, ptrNo bifurc_temp, Point A, Point B, Point C, 
                                int M, double *melhorCusto, bool *encontrou_valido,
                                double Qterm, double gamma, double mu, double epsilon, 
                                double R, Segment *segs_atuais, int total_segs);

// Seleciona o melhor segmento candidato minimizando o Volume Intravascular
ptrNo seleciona_melhor_cco1(Point novo, ptrNo raiz,
                            double R, double epsilon,
                            int M, double Qterm, double gamma, double mu,
                            int *rejeitados, Point *melhor_X_global);

#endif
