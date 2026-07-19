# Simulação de Difusão de Calor em Domínio 2D — MPI + OpenMP

Implementação híbrida da equação do calor em duas dimensões, desenvolvida para
a disciplina de IPPD.

## Execução

```sh
# compilar
mpicc -O2 -fopenmp -Wall -o heat2d heat2d.c -lm      # ou: make run

# executar: 4 processos MPI x 2 threads OpenMP cada
OMP_NUM_THREADS=2 mpirun -np 4 ./heat2d entrada.txt saida

# referência sequencial (para calcular o speedup)
OMP_NUM_THREADS=1 mpirun -np 1 ./heat2d entrada.txt saida_seq
```

Argumentos: `./heat2d [arquivo_de_entrada] [diretorio_de_saida]`
(padrões: `entrada.txt` e `saida`).

Ao final o programa imprime o tempo total (`MPI_Wtime`, máximo entre os
processos) e a taxa de células atualizadas por segundo — use esses números
para as medições de speedup/eficiência variando `-np` e `OMP_NUM_THREADS`.

## Visualização

```sh
python visualizar.py saida          # gera PNGs em saida/imagens/
python visualizar.py saida --gif    # também gera saida/animacao.gif
```

## Integrantes

Fabricio Bartz, Leonardo Braga e Victor Reis.
