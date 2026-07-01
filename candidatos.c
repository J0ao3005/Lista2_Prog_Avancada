#include "candidatos.h"

Point gera_ponto_aleatorio(double R) {
    Point p;
    do {
        p.x = ((double)rand() / RAND_MAX) * 2.0 * R - R;
        p.y = ((double)rand() / RAND_MAX) * 2.0 * R - R;
    } while (!dentro_dominio(p, R));
    return p;
}

void coleta_filhos(ptrNo raiz, ptrNo *lista, int *total) {
    if (raiz == NULL) return;
    if (raiz->pai != NULL) lista[(*total)++] = raiz;
    coleta_filhos(raiz->esq, lista, total);
    coleta_filhos(raiz->dir, lista, total);
}

void coleta_segmentos(ptrNo raiz, Segment *lista, int *total) {
    if (raiz == NULL) return;
    if (raiz->pai != NULL) {
        lista[*total].a = raiz->pai->p;
        lista[*total].b = raiz->p;
        (*total)++;
    }
    coleta_segmentos(raiz->esq, lista, total);
    coleta_segmentos(raiz->dir, lista, total);
}

// Parte F: Coordenadas Baricêntricas
Point pontoBaricentrico(Point A, Point B, Point C, double alpha, double beta, double lambda) {
    Point X;
    X.x = alpha * A.x + beta * B.x + lambda * C.x;
    X.y = alpha * A.y + beta * B.y + lambda * C.y;
    return X;
}

int pontoDentroTriangulo(Point A, Point B, Point C, Point P) {
    double o1 = orient(A, B, P);
    double o2 = orient(B, C, P);
    double o3 = orient(C, A, P);
    bool has_neg = (o1 < 0) || (o2 < 0) || (o3 < 0);
    bool has_pos = (o1 > 0) || (o2 > 0) || (o3 > 0);
    return !(has_neg && has_pos);
}

// Parte F e G: Busca na grade e avaliação de volume
Point otimizaBifurcacaoPorGrade(ptrNo raiz, ptrNo bifurc_temp, Point A, Point B, Point C, 
                                int M, double *melhorCusto, bool *encontrou_valido,
                                double Qterm, double gamma, double mu, double epsilon, 
                                double R, Segment *segs_atuais, int total_segs) {
    
    Point melhor_X = A;
    *encontrou_valido = false;

    // Discretização baseada na recomendação do roteiro
    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= M - i; j++) {
            
            double alpha  = (double)i / M;
            double beta   = (double)j / M;
            double lambda = 1.0 - alpha - beta;

            Point X = pontoBaricentrico(A, B, C, alpha, beta, lambda);

            // Garante que o ponto X da grade não vaze do domínio circular
            if (!dentro_dominio(X, R)) continue;

            // 1. Altera temporariamente a posição da bifurcação
            bifurc_temp->p = X;

            // 2. Prepara novos segmentos temporários para checagem de colisão
            Segment sAX = {A, X};
            Segment sXB = {X, B};
            Segment sXC = {X, C};

            bool colisao = false;
            for (int k = 0; k < total_segs; k++) {
                // Pula o segmento "pai -> filho" original, pois ele foi substituído
                if (fabs(segs_atuais[k].a.x - A.x) < EPSILON && fabs(segs_atuais[k].a.y - A.y) < EPSILON &&
                    fabs(segs_atuais[k].b.x - B.x) < EPSILON && fabs(segs_atuais[k].b.y - B.y) < EPSILON) {
                    continue;
                }

                // 6. Verifique se há interseção geométrica
                if (intersecta(sAX, segs_atuais[k]) || 
                    intersecta(sXB, segs_atuais[k]) || 
                    intersecta(sXC, segs_atuais[k])) {
                    colisao = true; 
                    break;
                }
                
                // Checagem de tolerância (epsilon)
                if (dist_ponto_segmento(X, segs_atuais[k]) < epsilon ||
                    dist_ponto_segmento(C, segs_atuais[k]) < epsilon) {
                    colisao = true; 
                    break;
                }
            }

            if (colisao) continue;

            // 2, 3 e 4. Atualiza comprimentos, fluxos e raios de TODA a árvore temporária
            atualizaGeometriaFisica(raiz, Qterm, gamma, mu);

            // 5. Calcula o volume total (Função Custo)
            double volume_atual = funcaoCustoVolume(raiz);

            // 7. Guarde o melhor valor de J
            if (volume_atual < *melhorCusto) {
                *melhorCusto = volume_atual;
                melhor_X = X;
                *encontrou_valido = true;
            }
        }
    }
    return melhor_X;
}

// Loop principal de seleção que itera sobre a árvore
ptrNo seleciona_melhor_cco1(Point novo, ptrNo raiz,
                            double R, double epsilon,
                            int M, double Qterm, double gamma, double mu,
                            int *rejeitados, Point *melhor_X_global) {

    Segment segs_atuais[MAX_NOS];
    int total_segs = 0;
    coleta_segmentos(raiz, segs_atuais, &total_segs);

    ptrNo filhos[MAX_NOS];
    int total_filhos = 0;
    coleta_filhos(raiz, filhos, &total_filhos);

    ptrNo melhor_filho = NULL;
    double menor_custo_global = DBL_MAX;
    bool rejeitado_nesta_rodada = true;

    for (int i = 0; i < total_filhos; i++) {
        ptrNo filho = filhos[i];
        ptrNo pai   = filho->pai;

        Point A = pai->p;
        Point B = filho->p;
        Point C = novo;

        // --- INÍCIO DA MODIFICAÇÃO TEMPORÁRIA ---
        No bifurc_temp = {0};
        No novo_no_temp = {0};

        bifurc_temp.pai = pai;
        bifurc_temp.esq = filho;
        bifurc_temp.dir = &novo_no_temp;

        novo_no_temp.pai = &bifurc_temp;
        novo_no_temp.p   = C;

        if (pai->esq == filho) pai->esq = &bifurc_temp;
        else pai->dir = &bifurc_temp;
        filho->pai = &bifurc_temp;
        // ----------------------------------------

        bool encontrou_valido = false;
        double melhor_custo_local = DBL_MAX;

        Point melhor_X_local = otimizaBifurcacaoPorGrade(
            raiz, &bifurc_temp, A, B, C, M, 
            &melhor_custo_local, &encontrou_valido, 
            Qterm, gamma, mu, epsilon, R, segs_atuais, total_segs
        );

        // --- DESFAZ A MODIFICAÇÃO TEMPORÁRIA ---
        if (pai->esq == &bifurc_temp) pai->esq = filho;
        else pai->dir = filho;
        filho->pai = pai;
        // ---------------------------------------

        if (encontrou_valido) {
            rejeitado_nesta_rodada = false;
            if (melhor_custo_local < menor_custo_global) {
                menor_custo_global = melhor_custo_local;
                melhor_filho = filho;
                *melhor_X_global = melhor_X_local;
            }
        }
    }

    // Restaura a física para o estado puro antes de retornar
    atualizaGeometriaFisica(raiz, Qterm, gamma, mu);

    if (rejeitado_nesta_rodada) {
        (*rejeitados)++;
    }

    return melhor_filho;
}
