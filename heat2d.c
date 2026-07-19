#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#ifdef _WIN32
#include <direct.h>
#define CRIAR_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define CRIAR_DIR(p) mkdir((p), 0755)
#endif

#define MAX_FONTES 64
#define LINHA_MAX 512

typedef struct
{
    int nx, ny;
    double alpha;
    double dx, dy, dt;
    int npassos;
    int salvar_a_cada;
    double t_topo, t_base;
    double t_esq, t_dir;
    double t_inicial;
    int nfontes;
} Parametros;

typedef struct
{
    int i0, j0, i1, j1;
    double temp;
} Fonte;

static int proxima_linha(FILE *f, char *buf)
{
    while (fgets(buf, LINHA_MAX, f))
    {
        char *p = buf;
        while (*p == ' ' || *p == '\t')
            p++;
        if (*p != '\0' && *p != '\n' && *p != '\r' && *p != '#')
            return 1;
    }
    return 0;
}

static int ler_entrada(const char *caminho, Parametros *par, Fonte *fontes)
{
    char buf[LINHA_MAX];
    FILE *f = fopen(caminho, "r");
    if (!f)
    {
        fprintf(stderr, "[rank 0] ERRO: Não foi possível abrir '%s'.\n", caminho);
        return 0;
    }

    if (!proxima_linha(f, buf) || sscanf(buf, "%d %d", &par->nx, &par->ny) != 2)
        goto erro;
    if (!proxima_linha(f, buf) || sscanf(buf, "%lf", &par->alpha) != 1)
        goto erro;
    if (!proxima_linha(f, buf) || sscanf(buf, "%lf %lf %lf",
                                         &par->dx, &par->dy, &par->dt) != 3)
        goto erro;
    if (!proxima_linha(f, buf) || sscanf(buf, "%d %d",
                                         &par->npassos, &par->salvar_a_cada) != 2)
        goto erro;
    if (!proxima_linha(f, buf) || sscanf(buf, "%lf %lf %lf %lf", &par->t_topo,
                                         &par->t_base, &par->t_esq, &par->t_dir) != 4)
        goto erro;
    if (!proxima_linha(f, buf) || sscanf(buf, "%lf", &par->t_inicial) != 1)
        goto erro;
    if (!proxima_linha(f, buf) || sscanf(buf, "%d", &par->nfontes) != 1)
        goto erro;

    if (par->nfontes < 0 || par->nfontes > MAX_FONTES)
    {
        fprintf(stderr, "[rank 0] ERRO: Número de fontes inválido (%d).\n", par->nfontes);
        fclose(f);
        return 0;
    }
    for (int k = 0; k < par->nfontes; k++)
    {
        if (!proxima_linha(f, buf) ||
            sscanf(buf, "%d %d %d %d %lf", &fontes[k].i0, &fontes[k].j0,
                   &fontes[k].i1, &fontes[k].j1, &fontes[k].temp) != 5)
            goto erro;
    }
    fclose(f);

    if (par->nx < 3 || par->ny < 3 || par->npassos < 1 || par->salvar_a_cada < 1)
    {
        fprintf(stderr, "[rank 0] ERRO: Parâmetros inválidos (nx/ny >= 3, passos >= 1).\n");
        return 0;
    }
    return 1;

erro:
    fprintf(stderr, "[rank 0] ERRO: Formato inválido no arquivo de entrada.\n");
    fclose(f);
    return 0;
}

static void montar_grade_inicial(const Parametros *par, const Fonte *fontes,
                                 double *grade)
{
    int nx = par->nx, ny = par->ny;

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            grade[i * ny + j] = par->t_inicial;

    for (int j = 0; j < ny; j++)
    {
        grade[0 * ny + j] = par->t_topo;
        grade[(nx - 1) * ny + j] = par->t_base;
    }
    for (int i = 0; i < nx; i++)
    {
        grade[i * ny + 0] = par->t_esq;
        grade[i * ny + (ny - 1)] = par->t_dir;
    }

    for (int k = 0; k < par->nfontes; k++)
    {
        for (int i = fontes[k].i0; i <= fontes[k].i1; i++)
        {
            if (i < 0 || i >= nx)
                continue;
            for (int j = fontes[k].j0; j <= fontes[k].j1; j++)
            {
                if (j < 0 || j >= ny)
                    continue;
                grade[i * ny + j] = fontes[k].temp;
            }
        }
    }
}

static void aplicar_celulas_fixas(const Parametros *par, const Fonte *fontes,
                                  double *T, int nlocal, int i_ini)
{
    int ny = par->ny;

    if (i_ini == 0)
        for (int j = 0; j < ny; j++)
            T[1 * ny + j] = par->t_topo;
    if (i_ini + nlocal == par->nx)
        for (int j = 0; j < ny; j++)
            T[nlocal * ny + j] = par->t_base;

    for (int k = 0; k < par->nfontes; k++)
    {
        int gi0 = fontes[k].i0 > i_ini ? fontes[k].i0 : i_ini;
        int gi1 = fontes[k].i1 < i_ini + nlocal - 1 ? fontes[k].i1 : i_ini + nlocal - 1;
        for (int gi = gi0; gi <= gi1; gi++)
        {
            int l = gi - i_ini + 1;
            for (int j = fontes[k].j0; j <= fontes[k].j1; j++)
            {
                if (j < 0 || j >= ny)
                    continue;
                T[l * ny + j] = fontes[k].temp;
            }
        }
    }
}

static void salvar_instante(const char *dir, int passo, double tempo,
                            const Parametros *par, const double *grade)
{
    char caminho[512];
    snprintf(caminho, sizeof caminho, "%s/passo_%06d.txt", dir, passo);

    FILE *f = fopen(caminho, "w");
    if (!f)
    {
        fprintf(stderr, "[rank 0] AVISO: Não foi possível gravar '%s'.\n", caminho);
        return;
    }
    fprintf(f, "# passo=%d tempo=%.6f nx=%d ny=%d\n", passo, tempo, par->nx, par->ny);
    for (int i = 0; i < par->nx; i++)
    {
        for (int j = 0; j < par->ny; j++)
            fprintf(f, "%.4f%c", grade[i * par->ny + j],
                    j == par->ny - 1 ? '\n' : ' ');
    }
    fclose(f);
}

int main(int argc, char **argv)
{
    int rank, nprocs, provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const char *arq_entrada = (argc > 1) ? argv[1] : "entrada.txt";
    const char *dir_saida = (argc > 2) ? argv[2] : "saida";

    Parametros par;
    Fonte fontes[MAX_FONTES];
    int ok = 1;

    if (rank == 0)
    {
        ok = ler_entrada(arq_entrada, &par, fontes);
        if (ok)
        {
            CRIAR_DIR(dir_saida);
            double fator = par.alpha * par.dt *
                           (1.0 / (par.dx * par.dx) + 1.0 / (par.dy * par.dy));
            if (fator > 0.5)
                fprintf(stderr, "AVISO: dt viola o critério de estabilidade "
                                "(alpha*dt*(1/dx^2+1/dy^2) = %.3f > 0.5); "
                                "a solução pode divergir.\n\n",
                        fator);
            printf("\nDifusão de Calor 2D: Grade %dx%d, %d Passos, "
                   "%d Processo(s) MPI x %d Thread(s) OpenMP\n\n",
                   par.nx, par.ny, par.npassos, nprocs, omp_get_max_threads());
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!ok)
    {
        MPI_Finalize();
        return 1;
    }
    MPI_Bcast(&par, sizeof(Parametros), MPI_BYTE, 0, MPI_COMM_WORLD);
    if (par.nfontes > 0)
        MPI_Bcast(fontes, par.nfontes * (int)sizeof(Fonte), MPI_BYTE, 0, MPI_COMM_WORLD);

    int nx = par.nx, ny = par.ny;

    int *contagens = malloc(nprocs * sizeof(int));
    int *deslocs = malloc(nprocs * sizeof(int));
    int base = nx / nprocs, resto = nx % nprocs, acum = 0;
    for (int r = 0; r < nprocs; r++)
    {
        int linhas = base + (r < resto ? 1 : 0);
        contagens[r] = linhas * ny;
        deslocs[r] = acum * ny;
        acum += linhas;
    }
    int nlocal = contagens[rank] / ny;
    int i_ini = deslocs[rank] / ny;

    if (nlocal == 0 && rank == 0)
        fprintf(stderr, "AVISO: Há mais processos que linhas na grade.\n");

    double *T = calloc((size_t)(nlocal + 2) * ny, sizeof(double));
    double *Tnew = calloc((size_t)(nlocal + 2) * ny, sizeof(double));
    double *grade_global = NULL;
    if (rank == 0)
        grade_global = malloc((size_t)nx * ny * sizeof(double));

    if (rank == 0)
        montar_grade_inicial(&par, fontes, grade_global);
    MPI_Scatterv(grade_global, contagens, deslocs, MPI_DOUBLE,
                 &T[ny], nlocal * ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    memcpy(Tnew, T, (size_t)(nlocal + 2) * ny * sizeof(double));

    int viz_cima = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
    int viz_baixo = (rank < nprocs - 1) ? rank + 1 : MPI_PROC_NULL;

    const double cx = par.alpha * par.dt / (par.dx * par.dx);
    const double cy = par.alpha * par.dt / (par.dy * par.dy);

    int l_ini = (i_ini == 0) ? 2 : 1;
    int l_fim = (i_ini + nlocal == nx) ? nlocal - 1 : nlocal;

    if (rank == 0)
        salvar_instante(dir_saida, 0, 0.0, &par, grade_global);

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    for (int passo = 1; passo <= par.npassos; passo++)
    {

        MPI_Sendrecv(&T[1 * ny], ny, MPI_DOUBLE, viz_cima, 0,
                     &T[(nlocal + 1) * ny], ny, MPI_DOUBLE, viz_baixo, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&T[nlocal * ny], ny, MPI_DOUBLE, viz_baixo, 1,
                     &T[0 * ny], ny, MPI_DOUBLE, viz_cima, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#pragma omp parallel for schedule(static)
        for (int l = l_ini; l <= l_fim; l++)
        {
            for (int j = 1; j < ny - 1; j++)
            {
                double c = T[l * ny + j];
                Tnew[l * ny + j] = c + cx * (T[(l + 1) * ny + j] - 2.0 * c + T[(l - 1) * ny + j]) + cy * (T[l * ny + (j + 1)] - 2.0 * c + T[l * ny + (j - 1)]);
            }
        }

        aplicar_celulas_fixas(&par, fontes, Tnew, nlocal, i_ini);

        // troca os buffers
        double *tmp = T;
        T = Tnew;
        Tnew = tmp;

        if (passo % par.salvar_a_cada == 0 || passo == par.npassos)
        {
            MPI_Gatherv(&T[ny], nlocal * ny, MPI_DOUBLE,
                        grade_global, contagens, deslocs, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
            if (rank == 0)
            {
                salvar_instante(dir_saida, passo, passo * par.dt, &par, grade_global);
                printf("passo %6d / %d  (t = %.4f)  ->  %s/passo_%06d.txt\n",
                       passo, par.npassos, passo * par.dt, dir_saida, passo);
            }
        }
    }

    double t1 = MPI_Wtime();
    double t_local = t1 - t0, t_max;
    MPI_Reduce(&t_local, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("\nTempo Total da Simulação: %.3f s "
               "(%d Processo(s) MPI x %d Thread(s) OpenMP)\n",
               t_max, nprocs, omp_get_max_threads());
        printf("Células atualizadas por Segundo: %.3e\n",
               (double)(nx - 2) * (ny - 2) * par.npassos / t_max);
    }

    free(T);
    free(Tnew);
    free(grade_global);
    free(contagens);
    free(deslocs);
    MPI_Finalize();
    return 0;
}
