# article06

Código de apoio ao artigo **Acoplamento Térmico–Fadiga em Interconexões BGA de
Processadores de Inteligência Artificial Orbitais: Estrutura Metodológica,
Contratividade do Modelo Térmico e Estimação de Vida por Observador**.

> **Aviso.** Os parâmetros do modelo térmico são **premissas ilustrativas**,
> declaradas no artigo (Seção de premissas). Não há dados experimentais, e os
> resultados **não são previsões de vida em missão**. O código demonstra a
> metodologia e a sensibilidade relativa.

## Estrutura

```
article06/
├── README.md
├── .gitignore
├── main.tex          # artigo
├── run_all.py        # executa os scripts e grava results/
├── script01.py       # modelo de referência (beta = 0): cenários, fase, Monte Carlo, EKF
├── script02.py       # invariância e contração (modelo de referência e estendido)
├── script03.py       # modelo estendido: modos de atitude e beta anual
└── results/          # saídas (JSON e manifest.json)
```

## Como executar

Requer apenas a biblioteca padrão do Python (3.8 ou superior; registre aqui a
versão usada: `python --version`).

```bash
python run_all.py                      # executa os 3 scripts
python run_all.py script03.py          # executa só um
python run_all.py --out outra_pasta    # muda a pasta de saída
```

Cada script define uma variável `result`. O `run_all.py` a lê e grava
`results/<script>_out.json`. Também grava `results/manifest.json` com o SHA-256
de cada script, a versão do Python e o tempo de execução. O código de saída é 0
se tudo rodou, e 1 se algum script falhou.

Rodando um script isoladamente (`python script03.py`), `script02.py` e
`script03.py` imprimem o JSON no terminal. O `script01.py` pode não imprimir
nada, e então só o `run_all.py` grava a saída dele.

## O que cada script gera

| Script | Conteúdo | Seções do artigo |
|---|---|---|
| `script01.py` | Órbita, cenários LEO/SSO idealizada, varredura de fase da carga, Monte Carlo (N=150, semente 2026), decomposição do viés do EKF, sensibilidade a erro de medição, reprodução algébrica do manuscrito original | Resultados do modelo de referência |
| `script02.py` | Condição de invariância, T_inv e constante de contração (λ) para o modelo de referência e para os modos zenith/side | Proposições 1 e 2 |
| `script03.py` | β ao longo do ano (SSO 6h, SSO 10h30, i = 51,6°), dias com eclipse, dano por órbita em função de β e vida anual por atitude × política de carga | Resultados do modelo estendido |

## Valores de conferência

Após rodar, os valores abaixo devem aparecer nos JSONs (tolerância de arredondamento):

| Arquivo | Campo | Esperado |
|---|---|---|
| `script01_out.json` | `orbita.T_min` | 94,6163 |
| `script01_out.json` | `cenarios.LEO_eclipse.vida_anos` | 9,1818 |
| `script01_out.json` | `ekf.C_ruido_e_erro_modelo.razao_media` | 0,3347 |
| `script02_out.json` | `referencia_beta0.lambda_W_K` | 0,4717 |
| `script02_out.json` | `estendido_zenith.fator_por_orbita` | 0,7481 |
| `script03_out.json` | `orbitas.SSO_10h30.vida_anos["zenith/mid"]` | 1,9319 |
| `script03_out.json` | `orbitas.i51.6_precessando.beta_max` | 74,816 |

Valores que divirjam além do arredondamento indicam alteração no código ou nos
parâmetros.

## Premissas

Todas as premissas estão na tabela de premissas do artigo. Em resumo:
- **Sem fonte** (ilustrativos): `Cc`, `Cs`, `G`, `eps`, `Ar`, `Poth`, potências, ruído e quantização dos sensores.
- **Do manuscrito original**: `F`, `L_DNP`, `h`, Δα.
- **Valores típicos sem fonte verificada**: `S0`, `Q_IR`, `ALPHA` (absortividade solar) e o fator de visada lateral do `script03.py`.

Troque pelos dados do seu hardware editando os dicionários no topo de cada script
e rode `python run_all.py` de novo.

## Limitações

- Modelo térmico de dois nós (die e PCB). Sem gradientes locais, albedo, penumbra, arrasto nem variação sazonal de S0.
- A atitude entra **somente como modo de apontamento** (condição de contorno térmica). Não há dinâmica de corpo rígido, inércia, torques de perturbação nem controle de atitude.
- Efemérides solares simplificadas; β calculado com erro da ordem de décimos de grau.
- Dano anual calculado em grade de 15 valores de β com interpolação linear. Ciclos lentos anuais não são contados.
- A relação de Engelmaier foi calibrada para outra faixa de temperatura. As temperaturas simuladas estão provavelmente fora do domínio.
- Sem prova de estabilidade do laço fechado com EKF.

## Como citar

Informe o **hash do commit** ou a **tag** usada no envio do artigo.
`[PENDENTE: preencher após criar a tag]`

## Licença

`[PENDENTE: escolher uma licença, por exemplo MIT]`
