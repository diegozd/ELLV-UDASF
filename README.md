# ELLV - UDASF
Algoritmo para cálculo de Equilíbrio-Líquido-Líquido-Vapor (ELLV) aplicado ao processo de desasfaltação (UDASF). Trabalho da disciplina de Termodinâmica de Soluções - COQ712 - PEQ/COPPE/UFRJ.

Aluno: Diego Telles Fernandes
PRofessores: Frederico W. Tavares & Iuri S. V. Segtovich

Este algoritmo calcula o rendimento de ODES e RASF a partir de correntes de Resíduo de Vácuo que foram Caracterizadas por Reconstrução Molecular. 

## Executando:

Para executar o código você deve criar uma pasta dentro do diretório **run** com o nome do RV que deseja avaliar. Esta pasta deve conter 3 arquivos. O arquivo que carrega a representação mais detalhada das moléculas da corrente de RV (como grafos) com nome de **Ligas.csv**, o arquivo que contem a representação simplificada das moléculas com nome **MSOL_expandida.csv** e o arquivo que contém a composição de cada molécula com nome de **Composicao.csv**. Veja exemplos na pasta **run** e consulte o relatório do trabalho na pasta **doc** para mais explicações.

Após ajeitada a estrutura de dados, o binário de nome UDASF contido na pasta **run** deve ser executado via terminal utilizando-se as seguintes flags: 

* `-T`      
Temperatura da extração em K.

* `-TC`      
Temperatura da extração em °C.

* `-P`      
Pressão da extração em kgf/cm²man.

* `-RSO`      
Razão Solvente Óleo (RSO) mássica.

* `-SOLV`      
Solvente utilizado para extração, podendo ser "C3" para propano ou "GLP". Caso deseje mudar alguma propriedade do solvente,  ou sua composição, alterar os arquivos das pastas relacionadas a cada solvente no diretório **run**.

* `-RV`      
Informar nome da pastas dentro do diretório **run** que contém os arquivos de Reconstrução Molecular que deseja simular.

#### Exemplo:
Na pasta raiz do projeto, execute:
``` ./run/UDASF -TC 60 -P 35 -RSO 6 -SOLV C3 -RV RV-A```

## Interpretando resultados:
Ao executar o binário o algoritmo criará uma corrente de RV, uma corrente de Solvente e uma corrente da mistura das duas primeiras. Posteriormente as fases da corrente de mistura são separadas e o solvente é removido de cada uma delas. Todas elas são avaliadas quanto a estabilidade, e possível separação de fases.  

O algoritmo mostrará estes resultados para cada corrente no terminal, indicando: pressão e temperatura; número de moléculas; número de fases; estado físico de cada fase (0-líquido; 1-vapor); fator de compressibilidade de cada fase; fração de cada fase; composição de cada fase; vazões mássicas e molares, a massa molar e a massa específica de cada fase, e da mistura global de todas as fases.

Depois de mostrar estas informações para cada corrente gerada, os resultados de rendimentos da desasfaltação é exibido, sendo estes os rendimentos de ODES e RASF em massa. A seção que indica os rendimentos é salva num arquivo dentro da própria pasta do RV. Em seguida será exibido um exemplo da resposta ao terminal.


<div style="overflow:scroll; width:600px; height:300px;">

~~~
ELLV para moleculas de Reconstrucao molecular
Diego Telles Fernandes
PEQ - COPPE - UFRJ
Trabalho da disciplina Termodinamica de Solucoes (COQ712) 
Compilado em Nov 10 2021 as 09:52:12

Wed Nov 10 10:24:52 2021



####    CORRENTE RV-A    ####

    RESULTADOS ELLV:    
T:      333.1 K
P:      35.3 bar
Nmol:   100
Nfases: 1
EstFis (0->Liquido ; 1->Vapor):
    0 
Z :
    1.0476 
beta:
    1.0000 
y:
00: 0.1527 
01: 0.0102 
02: 0.0001 
03: 0.1975 
04: 0.0712 
05: 0.0052 
06: 0.0019 
07: 0.0190 
08: 0.0001 
09: 0.0002 
10: 0.0003 
11: 0.0003 
12: 0.0288 
13: 0.0001 
14: 0.0004 
15: 0.0011 
16: 0.0022 
17: 0.0309 
18: 0.0024 
19: 0.0004 
20: 0.0177 
21: 0.0093 
22: 0.0060 
23: 0.0003 
24: 0.0042 
25: 0.0029 
26: 0.0005 
27: 0.0014 
28: 0.0008 
29: 0.0000 
30: 0.0112 
31: 0.0231 
32: 0.0118 
33: 0.0063 
34: 0.0000 
35: 0.0032 
36: 0.0131 
37: 0.0009 
38: 0.0162 
39: 0.0048 
40: 0.0000 
41: 0.0001 
42: 0.0095 
43: 0.0017 
44: 0.0019 
45: 0.0034 
46: 0.0011 
47: 0.0002 
48: 0.0040 
49: 0.0000 
50: 0.0046 
51: 0.0296 
52: 0.0056 
53: 0.0005 
54: 0.0003 
55: 0.0112 
56: 0.0034 
57: 0.0076 
58: 0.0002 
59: 0.0020 
60: 0.0015 
61: 0.0107 
62: 0.0010 
63: 0.0146 
64: 0.0003 
65: 0.0048 
66: 0.0027 
67: 0.0023 
68: 0.0002 
69: 0.0011 
70: 0.0377 
71: 0.0004 
72: 0.0000 
73: 0.0004 
74: 0.0062 
75: 0.0028 
76: 0.0018 
77: 0.0003 
78: 0.0085 
79: 0.0020 
80: 0.0037 
81: 0.0014 
82: 0.0001 
83: 0.0054 
84: 0.0076 
85: 0.0141 
86: 0.0004 
87: 0.0100 
88: 0.0010 
89: 0.0441 
90: 0.0171 
91: 0.0001 
92: 0.0006 
93: 0.0001 
94: 0.0209 
95: 0.0073 
96: 0.0009 
97: 0.0043 
98: 0.0028 
99: 0.0135 
Sum y:
1.0000 

    PROPRIEDADES DAS FASES (j) E DA MISTURA (mix):
j: 0 MMw:786.4051        roh:957.6035    VazN:1.0000     VazM:786.4051
mix: MMw:786.4051        roh:957.6035    VazN:1.0000     VazM:786.4051
 Objeto RV-A construido com sucesso



####    CORRENTE C3    ####

    RESULTADOS ELLV:    
T:      333.1 K
P:      35.3 bar
Nmol:   1
Nfases: 1
EstFis (0->Liquido ; 1->Vapor):
    0 
Z :
    0.1279 
beta:
    1.0000 
y:
00: 1.0000 
Sum y:
1.0000 

    PROPRIEDADES DAS FASES (j) E DA MISTURA (mix):
j: 0 MMw:44.0000         roh:438.8886    VazN:71.4914    VazM:3145.6204
mix: MMw:44.0000         roh:438.8886    VazN:71.4914    VazM:3145.6204
 Objeto C3 construido com sucesso
 Real RSO SOLV:4



####    CORRENTE RVSOLV_MX    ####

    RESULTADOS ELLV:    
T:      333.1 K
P:      35.3 bar
Nmol:   101
Nfases: 2
EstFis (0->Liquido ; 1->Vapor):
    0           0 
Z :
    0.1321      0.3956 
beta:
    0.9813      0.0187 
y:
00: 0.9920      0.6829 
01: 0.0021      0.0046 
02: 0.0000      0.0051 
03: 0.0000      0.0000 
04: 0.0019      0.0483 
05: 0.0010      0.0021 
06: 0.0000      0.0032 
07: 0.0000      0.0014 
08: 0.0000      0.0136 
09: 0.0000      0.0000 
10: 0.0000      0.0001 
11: 0.0000      0.0002 
12: 0.0000      0.0002 
13: 0.0000      0.0192 
14: 0.0000      0.0000 
15: 0.0000      0.0002 
16: 0.0000      0.0007 
17: 0.0000      0.0011 
18: 0.0004      0.0007 
19: 0.0000      0.0011 
20: 0.0000      0.0001 
21: 0.0002      0.0009 
22: 0.0001      0.0038 
23: 0.0001      0.0013 
24: 0.0000      0.0001 
25: 0.0000      0.0006 
26: 0.0000      0.0003 
27: 0.0000      0.0001 
28: 0.0000      0.0003 
29: 0.0000      0.0006 
30: 0.0000      0.0000 
31: 0.0001      0.0051 
32: 0.0001      0.0102 
33: 0.0001      0.0018 
34: 0.0000      0.0037 
35: 0.0000      0.0000 
36: 0.0000      0.0019 
37: 0.0000      0.0083 
38: 0.0000      0.0001 
39: 0.0000      0.0106 
40: 0.0000      0.0017 
41: 0.0000      0.0000 
42: 0.0000      0.0000 
43: 0.0001      0.0042 
44: 0.0000      0.0003 
45: 0.0000      0.0006 
46: 0.0000      0.0015 
47: 0.0000      0.0008 
48: 0.0000      0.0000 
49: 0.0000      0.0029 
50: 0.0000      0.0000 
51: 0.0001      0.0007 
52: 0.0001      0.0155 
53: 0.0001      0.0014 
54: 0.0000      0.0002 
55: 0.0000      0.0000 
56: 0.0001      0.0019 
57: 0.0000      0.0013 
58: 0.0001      0.0015 
59: 0.0000      0.0000 
60: 0.0000      0.0009 
61: 0.0000      0.0010 
62: 0.0000      0.0075 
63: 0.0000      0.0006 
64: 0.0001      0.0080 
65: 0.0000      0.0001 
66: 0.0000      0.0031 
67: 0.0000      0.0019 
68: 0.0000      0.0004 
69: 0.0000      0.0001 
70: 0.0000      0.0005 
71: 0.0001      0.0240 
72: 0.0000      0.0003 
73: 0.0000      0.0000 
74: 0.0000      0.0003 
75: 0.0001      0.0008 
76: 0.0000      0.0013 
77: 0.0000      0.0011 
78: 0.0000      0.0002 
79: 0.0001      0.0016 
80: 0.0000      0.0012 
81: 0.0000      0.0027 
82: 0.0000      0.0002 
83: 0.0000      0.0000 
84: 0.0000      0.0025 
85: 0.0000      0.0044 
86: 0.0001      0.0072 
87: 0.0000      0.0003 
88: 0.0000      0.0072 
89: 0.0000      0.0004 
90: 0.0001      0.0283 
91: 0.0001      0.0052 
92: 0.0000      0.0001 
93: 0.0000      0.0004 
94: 0.0000      0.0000 
95: 0.0001      0.0087 
96: 0.0000      0.0054 
97: 0.0000      0.0001 
98: 0.0000      0.0013 
99: 0.0000      0.0004 
100: 0.0002     0.0007 
Sum y:
1.0000  1.0000 

    PROPRIEDADES DAS FASES (j) E DA MISTURA (mix):
j: 0 MMw:49.9881         roh:482.5769    VazN:71.1354    VazM:3555.9285
j: 1 MMw:277.3704        roh:894.5434    VazN:1.3559     VazM:376.0970
mix: MMw:54.2413         roh:504.8139    VazN:72.4914    VazM:3932.0255
 Objeto RVSOLV_MX construido com sucesso

 Objeto com vetor de correntes a partir de RVSOLV_MX foi construido com sucesso


####    CORRENTE ODES    ####

    RESULTADOS ELLV:    
T:      333.1 K
P:      35.3 bar
Nmol:   100
Nfases: 1
EstFis (0->Liquido ; 1->Vapor):
    0 
Z :
    1.0689 
beta:
    1.0000 
y:
00: 0.2570 
01: 0.0059 
02: 0.0001 
03: 0.2316 
04: 0.1198 
05: 0.0015 
06: 0.0000 
07: 0.0009 
08: 0.0001 
09: 0.0002 
10: 0.0000 
11: 0.0000 
12: 0.0050 
13: 0.0001 
14: 0.0002 
15: 0.0003 
16: 0.0012 
17: 0.0525 
18: 0.0014 
19: 0.0004 
20: 0.0288 
21: 0.0072 
22: 0.0075 
23: 0.0003 
24: 0.0058 
25: 0.0043 
26: 0.0006 
27: 0.0019 
28: 0.0000 
29: 0.0001 
30: 0.0075 
31: 0.0164 
32: 0.0163 
33: 0.0022 
34: 0.0000 
35: 0.0010 
36: 0.0031 
37: 0.0013 
38: 0.0031 
39: 0.0043 
40: 0.0000 
41: 0.0001 
42: 0.0066 
43: 0.0021 
44: 0.0018 
45: 0.0025 
46: 0.0000 
47: 0.0002 
48: 0.0001 
49: 0.0001 
50: 0.0064 
51: 0.0151 
52: 0.0065 
53: 0.0003 
54: 0.0005 
55: 0.0152 
56: 0.0028 
57: 0.0098 
58: 0.0003 
59: 0.0013 
60: 0.0002 
61: 0.0009 
62: 0.0001 
63: 0.0065 
64: 0.0002 
65: 0.0011 
66: 0.0001 
67: 0.0031 
68: 0.0000 
69: 0.0007 
70: 0.0090 
71: 0.0000 
72: 0.0000 
73: 0.0000 
74: 0.0091 
75: 0.0017 
76: 0.0006 
77: 0.0000 
78: 0.0112 
79: 0.0005 
80: 0.0000 
81: 0.0020 
82: 0.0001 
83: 0.0036 
84: 0.0029 
85: 0.0076 
86: 0.0000 
87: 0.0004 
88: 0.0008 
89: 0.0102 
90: 0.0177 
91: 0.0000 
92: 0.0000 
93: 0.0002 
94: 0.0159 
95: 0.0000 
96: 0.0012 
97: 0.0045 
98: 0.0039 
99: 0.0220 
Sum y:
1.0000 

    PROPRIEDADES DAS FASES (j) E DA MISTURA (mix):
j: 0 MMw:791.2798        roh:944.3594    VazN:0.5700     VazM:451.0506
mix: MMw:791.2798        roh:944.3594    VazN:0.5700     VazM:451.0506
VCRVSOLV.size():2


####    CORRENTE RASF    ####

    RESULTADOS ELLV:    
T:      333.1 K
P:      35.3 bar
Nmol:   100
Nfases: 1
EstFis (0->Liquido ; 1->Vapor):
    0 
Z :
    1.0198 
beta:
    1.0000 
y:
00: 0.0145 
01: 0.0160 
02: 0.0001 
03: 0.1523 
04: 0.0067 
05: 0.0101 
06: 0.0044 
07: 0.0429 
08: 0.0001 
09: 0.0002 
10: 0.0008 
11: 0.0006 
12: 0.0604 
13: 0.0000 
14: 0.0007 
15: 0.0022 
16: 0.0034 
17: 0.0022 
18: 0.0036 
19: 0.0005 
20: 0.0028 
21: 0.0121 
22: 0.0040 
23: 0.0003 
24: 0.0020 
25: 0.0009 
26: 0.0003 
27: 0.0009 
28: 0.0018 
29: 0.0000 
30: 0.0160 
31: 0.0321 
32: 0.0058 
33: 0.0118 
34: 0.0000 
35: 0.0061 
36: 0.0263 
37: 0.0004 
38: 0.0336 
39: 0.0055 
40: 0.0000 
41: 0.0000 
42: 0.0133 
43: 0.0010 
44: 0.0019 
45: 0.0046 
46: 0.0025 
47: 0.0001 
48: 0.0091 
49: 0.0000 
50: 0.0022 
51: 0.0488 
52: 0.0044 
53: 0.0006 
54: 0.0001 
55: 0.0060 
56: 0.0041 
57: 0.0046 
58: 0.0001 
59: 0.0030 
60: 0.0032 
61: 0.0237 
62: 0.0020 
63: 0.0253 
64: 0.0004 
65: 0.0098 
66: 0.0061 
67: 0.0012 
68: 0.0004 
69: 0.0016 
70: 0.0757 
71: 0.0008 
72: 0.0000 
73: 0.0009 
74: 0.0025 
75: 0.0043 
76: 0.0034 
77: 0.0007 
78: 0.0049 
79: 0.0039 
80: 0.0085 
81: 0.0006 
82: 0.0001 
83: 0.0077 
84: 0.0140 
85: 0.0226 
86: 0.0009 
87: 0.0228 
88: 0.0013 
89: 0.0891 
90: 0.0163 
91: 0.0003 
92: 0.0014 
93: 0.0001 
94: 0.0275 
95: 0.0170 
96: 0.0005 
97: 0.0041 
98: 0.0014 
99: 0.0023 
Sum y:
1.0000 

    PROPRIEDADES DAS FASES (j) E DA MISTURA (mix):
j: 0 MMw:779.9426        roh:975.7008    VazN:0.4300     VazM:335.3545
mix: MMw:779.9426        roh:975.7008    VazN:0.4300     VazM:335.3545

###############################################
####       RENDIMENTOS DESASFALTACAO       ####
###############################################

T: 333.15 K  =>  60.00 C
P: 35.34 bar =>  35.00 kgf/cm2g

OLEO:      RV-A
SOLVENTE:  C3
RSO:       4.00
rend ODES: 57.36 %
rend RASF: 42.64 %


Tempo de execucao: 37.34s

FIM!
~~~
</div>

## Instalação:

O único pré-requisito externo é a biblioteca Eigen, que já está inserida no projeto, então nada precisa ser instalado. 

Se você está em ambiente *linux* então basta clonar o projeto e executar o binário conforme exemplo anterior. Caso deseje recompilar o código, seu sistema operacional já deve possuir um compilador C++, então basta executar `make` na raiz do projeto. Estou utilizando o c++20 como compilador (`g++ -std=c++2a`), e caso deseje alterar a versão do C++, alterar alguma flag, ou adicionar alguma bibliotecas, faça isso no arquivo *makefile*. No entanto, lembre-se: **Duas classe deste código possuem restrições de propriedade intelectual e não estão disponíveis aqui no git**, caso queira recompilar, certifique-se de ter conseguido estes arquivos com os proprietários. 

Se está em ambiente *windows*, é possível adaptar o makefile para seu sistema operacional, mas pensaria melhor sobre utilizar algo menos proprietário para este fim. De qualquer modo, boa sorte. 

