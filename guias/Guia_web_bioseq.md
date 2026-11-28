# Guía paso a paso de la web bioseq

*Qué hace cada parte, para qué sirve y la biología que hay detrás.*

**[▶ Abrir la web](https://dcn-bioinformatics.github.io/BasicBioSet/web/)** · [Read this guide in English](Bioseq_web_guide.md)

## Índice

- [1. Qué es bioseq](#1-qué-es-bioseq)
- [2. Conceptos básicos de biología](#2-conceptos-básicos-de-biología)
  - [2.1 El ADN](#21-el-adn)
  - [2.2 El ARN](#22-el-arn)
  - [2.3 Genes, codones y proteínas](#23-genes-codones-y-proteínas)
  - [2.4 Los formatos de archivo: FASTA y FASTQ](#24-los-formatos-de-archivo-fasta-y-fastq)
- [3. Primeros pasos con la web](#3-primeros-pasos-con-la-web)
  - [3.1 Cómo abrir la web](#31-cómo-abrir-la-web)
  - [3.2 Aspecto general](#32-aspecto-general)
  - [3.3 Panel «Your sequences» (Tus secuencias)](#33-panel-your-sequences-tus-secuencias)
  - [3.4 Botones comunes en los resultados](#34-botones-comunes-en-los-resultados)
- [4. Las pestañas de análisis, una por una](#4-las-pestañas-de-análisis-una-por-una)
  - [Pestaña «Overview» (Resumen)](#pestaña-overview-resumen)
  - [Pestaña «Translation» (Traducción)](#pestaña-translation-traducción)
  - [Pestaña «ORFs» (Marcos abiertos de lectura)](#pestaña-orfs-marcos-abiertos-de-lectura)
  - [Pestaña «Codon usage» (Uso de codones)](#pestaña-codon-usage-uso-de-codones)
  - [Pestaña «Motif search» (Búsqueda de motivos)](#pestaña-motif-search-búsqueda-de-motivos)
  - [Pestaña «K-mers» (Palabras de longitud k)](#pestaña-k-mers-palabras-de-longitud-k)
  - [Pestaña «Composition» (Composición)](#pestaña-composition-composición)
  - [Pestaña «Tm & weight» (Temperatura de fusión y peso)](#pestaña-tm--weight-temperatura-de-fusión-y-peso)
  - [Pestaña «Read quality» (Calidad de las lecturas)](#pestaña-read-quality-calidad-de-las-lecturas)
  - [Pestaña «All sequences» (Todas las secuencias)](#pestaña-all-sequences-todas-las-secuencias)
- [5. Ejemplos de uso de principio a fin](#5-ejemplos-de-uso-de-principio-a-fin)
  - [5.1 Analizar un gen desconocido](#51-analizar-un-gen-desconocido)
  - [5.2 Preparar una clonación](#52-preparar-una-clonación)
  - [5.3 Revisar los datos de una secuenciación](#53-revisar-los-datos-de-una-secuenciación)
  - [5.4 Archivos de ejemplo incluidos](#54-archivos-de-ejemplo-incluidos)
- [6. Glosario](#6-glosario)

## 1. Qué es bioseq

**bioseq** es una caja de herramientas para analizar secuencias de ADN y ARN. La web permite pegar o subir una secuencia y obtener al instante su composición, su traducción a proteína, los genes posibles que contiene, los sitios de corte de enzimas, su temperatura de fusión y, si son lecturas de un secuenciador, la calidad de esas lecturas.

Todo el análisis se hace **dentro de tu navegador**: la web descarga una versión de Python (Pyodide) y ejecuta en tu propio ordenador el mismo código que usa la versión de terminal. Por eso:

- **Tus secuencias nunca salen de tu ordenador.** No se suben a ningún servidor (lo recuerda el pie de la página).
- La web y el menú de terminal (`python3 main.py`) dan **exactamente los mismos resultados**.
- La primera visita descarga unos 10 MB y tarda un poco; las siguientes son más rápidas.

La interfaz de la web está en inglés. En esta guía cada botón y pestaña aparece con su nombre original en inglés y su traducción entre paréntesis.

## 2. Conceptos básicos de biología

Esta sección explica lo mínimo necesario para entender lo que muestra la web. Si ya conoces el tema, puedes saltar a la sección 3.

### 2.1 El ADN

El **ADN** (ácido desoxirribonucleico) es la molécula que guarda la información genética de casi todos los seres vivos. Es una cadena larga formada por cuatro piezas llamadas **nucleótidos** o **bases**:

| Letra | Base | Se empareja con | Tipo |
|---|---|---|---|
| A | Adenina | T (2 puentes de hidrógeno) | Purina |
| T | Timina | A | Pirimidina |
| G | Guanina | C (3 puentes de hidrógeno) | Purina |
| C | Citosina | G | Pirimidina |

El ADN tiene **dos hebras** enrolladas en una doble hélice. Cada base de una hebra se empareja con su complementaria en la otra (A con T, G con C). Por eso, conociendo una hebra se conoce la otra. Las hebras tienen dirección: se leen de **5′ a 3′** (se pronuncia «cinco prima a tres prima»), y la hebra opuesta va en sentido contrario.

La longitud del ADN se mide en **pares de bases (bp)** o, para una sola hebra, en **nucleótidos (nt)**. Además de A, C, G y T, la web acepta la letra **N**, que significa «base desconocida» (por ejemplo, una posición que el secuenciador no pudo leer).

### 2.2 El ARN

El **ARN** es una copia de trabajo de un fragmento de ADN. Es parecido, pero usa **U (uracilo)** en lugar de T. El paso de ADN a ARN se llama **transcripción**. El ARN mensajero (**ARNm**) es el que lleva las instrucciones para fabricar una proteína.

La web acepta secuencias de ARN: si pegas una secuencia con U, la lee automáticamente como ADN (cambia U por T) y la marca como «RNA» en la lista de secuencias.

### 2.3 Genes, codones y proteínas

Las proteínas son cadenas de **aminoácidos** (hay 20 tipos). El ARNm se lee de tres en tres bases; cada grupo de tres es un **codón** y corresponde a un aminoácido. Esta correspondencia es el **código genético**:

- Hay 4 × 4 × 4 = **64 codones** posibles. 61 codifican aminoácidos y 3 son **codones de parada** (TAA, TAG, TGA), que indican «fin de la proteína».
- El codón **ATG** codifica metionina (M) y además marca el **inicio** de la mayoría de los genes.
- Varios codones pueden dar el mismo aminoácido (por ejemplo, la leucina tiene 6). Se dice que el código es «degenerado».
- Este paso de ARN a proteína se llama **traducción**. En la web los aminoácidos se escriben con su código de una letra (M, A, L…) y las paradas con un asterisco `*`.

Como los codones se leen de tres en tres, hay tres maneras de «cortar» una secuencia en codones según se empiece en la base 1, 2 o 3. Cada una es un **marco de lectura**. Si además se lee la hebra opuesta, hay **seis marcos** en total.

### 2.4 Los formatos de archivo: FASTA y FASTQ

**FASTA** es el formato más común para guardar secuencias. Cada secuencia empieza con una línea de cabecera que comienza por `>` seguida del nombre, y debajo va la secuencia:

```
>mi_gen descripción opcional
ATGGCTAGCAAAGGAGAAGAACTTTTC...
```

**FASTQ** es el formato que producen las máquinas de secuenciación. Cada **lectura** (fragmento de ADN leído por la máquina) ocupa 4 líneas: el nombre (empieza con `@`), la secuencia, una línea con `+` y una línea de **calidad** con un carácter por cada base. Ese carácter codifica cuán segura está la máquina de haber leído bien esa base (ver la pestaña «Read quality» en la sección 4).

## 3. Primeros pasos con la web

### 3.1 Cómo abrir la web

1. Abre una terminal en la carpeta del proyecto y ejecuta `python3 main.py --web`.
2. Abre en el navegador la dirección `http://127.0.0.1:8000/web/`.
3. Espera a que el indicador de estado, arriba a la derecha, deje de decir **Starting…** (Iniciando). Cuando el motor está listo, los botones se activan.

Si el puerto 8000 está ocupado, usa `python3 main.py --web --port 8080`. Para que otros ordenadores de tu red puedan entrar, añade `--host 0.0.0.0`. Si el proyecto está publicado en GitHub Pages, basta con abrir su dirección `.../web/`.

> **Importante:** Abrir directamente el archivo `web/index.html` con doble clic no funciona: el navegador bloquea la carga de los archivos de Python. Siempre hay que abrirla a través de un servidor (con `--web` o GitHub Pages).

### 3.2 Aspecto general

La pantalla tiene tres zonas:

- **Barra superior:** el nombre bioseq y el indicador de estado del motor de análisis.
- **Panel izquierdo «Your sequences» (Tus secuencias):** donde se introducen los datos y se elige qué secuencia analizar.
- **Panel derecho:** las pestañas de análisis, los controles de cada pestaña y los resultados.

### 3.3 Panel «Your sequences» (Tus secuencias)

Aquí se cargan los datos. Puedes introducirlos de tres formas:

| Elemento | Qué hace |
|---|---|
| Cuadro de texto | Pega una secuencia suelta (solo letras), un texto FASTA con una o varias secuencias, o un texto FASTQ. También puedes **arrastrar y soltar** un archivo encima. |
| **Analyze** (Analizar) | Lee el texto del cuadro y carga las secuencias. Los números y espacios de una secuencia suelta se ignoran, así que puedes pegar secuencias copiadas de GenBank con su numeración. |
| **Upload file** (Subir archivo) | Abre un archivo de tu ordenador. Acepta .fasta, .fa, .fna, .ffn, .fas, .seq, .fastq, .fq y .txt. |
| **Example** (Ejemplo) | Carga una secuencia de ejemplo para probar la web sin datos propios. |
| **Clear** (Limpiar) | Vacía el cuadro de texto. |
| Lista **Loaded** (Cargadas) | Muestra cada secuencia cargada con su nombre, longitud y % de GC. Haz clic en una para que las pestañas la analicen. Las de ARN llevan la etiqueta «RNA». |

Si el texto tiene un error (por ejemplo, letras que no son bases, o texto antes de la primera línea `>`), aparece un mensaje debajo de los botones explicando qué corregir. Solo se admiten A, C, G, T (o U) y N.

Con un archivo **FASTQ**, la web carga hasta 20 000 lecturas (para que el navegador no se ralentice) y abre directamente la pestaña «Read quality». Cada lectura aparece también en la lista, para poder analizarlas una a una con las demás pestañas.

### 3.4 Botones comunes en los resultados

- **Copy** (Copiar): copia una secuencia al portapapeles.
- **FASTA / CSV:** descargan el resultado. Los CSV se abren con Excel, LibreOffice o Google Sheets. Las descargas siempre contienen los datos completos aunque la pantalla muestre una versión recortada (las secuencias de más de 12 000 caracteres se muestran truncadas).
- **Mapas y gráficos:** al pasar el ratón por encima de una barra, un punto o una marca, aparece una etiqueta con el detalle.

## 4. Las pestañas de análisis, una por una

Las ocho primeras pestañas analizan **la secuencia seleccionada** en la lista de la izquierda. Las dos últimas («Read quality» y «All sequences») trabajan con **todas las secuencias cargadas** a la vez.

| Pestaña | En una frase |
|---|---|
| Overview (Resumen) | Longitud, cuántas A, C, G, T tiene, % de GC, hebra complementaria y ARNm. |
| Translation (Traducción) | Convierte el ADN en proteína en uno o en los seis marcos de lectura. |
| ORFs (Marcos abiertos de lectura) | Busca los posibles genes: tramos de ATG a un codón de parada. |
| Codon usage (Uso de codones) | Cuenta qué codones usan los genes encontrados para cada aminoácido. |
| Motif search (Búsqueda de motivos) | Encuentra una secuencia corta o el sitio de corte de una enzima. |
| K-mers | Cuenta las «palabras» de k letras más frecuentes. |
| Composition (Composición) | Gráficos de GC, sesgo GC y complejidad a lo largo de la secuencia; islas CpG. |
| Tm & weight (Tm y peso) | Temperatura de fusión de un cebador y peso molecular. |
| Read quality (Calidad de lecturas) | Calidad de las lecturas de un secuenciador (FASTQ) y su recorte. |
| All sequences (Todas las secuencias) | Tabla resumen de todas las secuencias cargadas. |

### Pestaña «Overview» (Resumen)

#### ¿Qué hace?

Es la ficha básica de la secuencia: cuántas bases tiene, de qué tipo, su contenido en GC, su hebra complementaria y su versión en ARN mensajero.

#### ¿Para qué sirve?

- Comprobar de un vistazo que la secuencia cargada es la correcta (longitud esperada, sin N inesperadas).
- Obtener la **hebra complementaria inversa**, que se necesita, por ejemplo, para diseñar el cebador «reverse» de una PCR o para buscar un gen que está en la otra hebra.
- Obtener el **ARNm** correspondiente.

#### Un poco de biología

El **contenido en GC** es el porcentaje de bases que son G o C. Es importante porque el par G–C está unido por 3 puentes de hidrógeno y el A–T solo por 2: cuanto más GC, más estable es el ADN y más temperatura hace falta para separar sus hebras. Cada especie tiene un GC típico (el ser humano ronda el 41 %, algunas bacterias superan el 70 %), así que también sirve para reconocer de dónde viene un fragmento.

La **complementaria inversa** es la otra hebra del ADN escrita también de 5′ a 3′: se invierte la secuencia y se cambia cada base por su pareja (A↔T, G↔C).

#### Cómo usarla, paso a paso

1. Carga una secuencia y selecciónala en la lista de la izquierda.
2. Abre la pestaña **Overview** (es la que aparece por defecto).

#### Cómo leer los resultados

- **Length:** longitud en pares de bases (bp).
- **GC content:** porcentaje de G + C.
- **A, C, G, T (y N si hay):** número de cada base y su porcentaje. La barra de colores muestra la proporción de cada una.
- **Reverse complement:** la hebra complementaria inversa, con botón para copiarla o descargarla en FASTA.
- **mRNA:** la secuencia transcrita (las T pasan a U).

### Pestaña «Translation» (Traducción)

#### ¿Qué hace?

Traduce la secuencia de ADN a proteína usando el código genético estándar, en el marco de lectura que elijas o en los seis a la vez.

#### ¿Para qué sirve?

- Saber qué proteína codifica un gen.
- Encontrar en qué marco de lectura está el gen: el marco correcto suele ser el que da una proteína larga sin asteriscos (paradas) en medio.
- Comprobar si una mutación cambia un aminoácido o introduce una parada prematura.

#### Un poco de biología

Recuerda (sección 2.3) que la secuencia se lee de tres en tres bases. Los marcos **+1, +2 y +3** empiezan a leer en la base 1, 2 y 3 de la secuencia tal como está escrita; los marcos **−1, −2 y −3** hacen lo mismo con la hebra complementaria inversa. Un gen real solo está en uno de ellos.

#### Cómo usarla, paso a paso

1. Selecciona la secuencia y abre **Translation**.
2. En **Reading frame** (Marco de lectura) elige **+1**, **+2**, **+3** o **All six** (los seis).
3. Opcional: marca **Stop at first stop codon** (Detenerse en el primer codón de parada) para que la proteína termine en la primera parada, como ocurre en la célula. Esta opción no está disponible al ver los seis marcos.

#### Cómo leer los resultados

- La proteína aparece con el código de una letra de los aminoácidos.
- Los codones de parada se muestran como un **asterisco rojo** (`*`).
- Los codones con alguna N, que no se pueden traducir, aparecen como **X**.
- En modo seis marcos se puede descargar todo en un único archivo FASTA.

### Pestaña «ORFs» (Marcos abiertos de lectura)

#### ¿Qué hace?

Busca automáticamente los **ORF** (del inglés «Open Reading Frame», marco abierto de lectura): tramos que empiezan en un codón ATG y terminan en un codón de parada, en los tres marcos de las dos hebras.

#### ¿Para qué sirve?

- Localizar posibles **genes** en una secuencia desconocida (por ejemplo, el genoma de un virus).
- Saber dónde empieza y termina cada gen, en qué hebra está y qué proteína produciría.
- Obtener las proteínas en FASTA para buscarlas después en bases de datos como BLAST o UniProt.

#### Un poco de biología

Un gen que codifica una proteína es, en lo esencial, un ORF: empieza en ATG y sigue sin interrupciones hasta una parada. Por azar aparecen muchos ORF muy cortos que no son genes reales; por eso se fija una longitud mínima. Cuanto más largo es un ORF, más probable es que sea un gen de verdad (en genomas de eucariotas, con intrones, los genes están partidos y esta búsqueda simple no los detecta enteros).

#### Cómo usarla, paso a paso

1. Selecciona la secuencia y abre **ORFs**.
2. Ajusta **Minimum ORF length** (Longitud mínima del ORF) en nucleótidos; por defecto 30 nt (10 aminoácidos). Para genomas reales conviene subirlo a 300 nt (100 aminoácidos) o más.

#### Cómo leer los resultados

- Un **mapa** de la secuencia con cada ORF dibujado como una flecha sobre la hebra + o − (pasa el ratón por encima para ver los detalles).
- Una **tabla** con: posición inicial y final (contando desde 1 sobre la hebra directa), hebra (+ o −), longitud en nucleótidos, longitud de la proteína en aminoácidos, su **peso molecular en kDa** y la secuencia de la proteína (clic para desplegarla si es larga).
- Botones para descargar la tabla en **CSV** y las proteínas en **FASTA**.
- Con muchos ORF, la tabla muestra los 500 primeros; el CSV los contiene todos.
- Si no aparece ninguno, la web lo indica y sugiere bajar la longitud mínima.

### Pestaña «Codon usage» (Uso de codones)

#### ¿Qué hace?

Cuenta cuántas veces aparece cada uno de los 64 codones dentro de los ORF encontrados (con la misma longitud mínima de la pestaña ORFs) y los agrupa por aminoácido.

#### ¿Para qué sirve?

- Ver qué codones «prefiere» un organismo para cada aminoácido.
- **Optimizar genes** para expresarlos en otro organismo (por ejemplo, producir una proteína humana en la bacteria «E. coli»): conviene usar los codones que el organismo de destino usa más.
- Detectar genes adquiridos de otra especie (transferencia horizontal), que suelen tener un uso de codones distinto al del resto del genoma.

#### Un poco de biología

Como varios codones codifican el mismo aminoácido (codones sinónimos), cada especie los usa con frecuencias distintas, en parte porque tiene más ARN de transferencia para unos que para otros. Esta preferencia se llama **sesgo de uso de codones**.

#### Cómo usarla, paso a paso

1. Selecciona la secuencia y abre **Codon usage**.
2. Ajusta la **longitud mínima del ORF** si hace falta (es el mismo valor que en la pestaña ORFs).

#### Cómo leer los resultados

- Una fila por aminoácido (nombre de tres letras y código de una letra) con el total de codones contados. Las paradas también se incluyen.
- Para cada codón: cuántas veces aparece y qué **porcentaje** representa entre los codones de ese aminoácido, con una barra. Los codones que no aparecen se muestran atenuados.
- Descarga en **CSV**.

### Pestaña «Motif search» (Búsqueda de motivos)

#### ¿Qué hace?

Busca todas las apariciones de una secuencia corta (un **motivo**) en las dos hebras. También reconoce los nombres de **enzimas de restricción** y de algunos motivos comunes, y en ese caso indica dónde corta la enzima.

#### ¿Para qué sirve?

- Planificar una **clonación**: saber si una enzima corta tu secuencia, cuántas veces y dónde.
- Encontrar señales reguladoras conocidas, como la **caja TATA** de los promotores.
- Buscar cualquier secuencia de interés (un cebador, una mutación, una etiqueta).

#### Un poco de biología

Las **enzimas de restricción** son proteínas de bacterias que cortan el ADN solo donde encuentran una secuencia concreta, su **sitio de reconocimiento** (por ejemplo, EcoRI reconoce GAATTC y corta entre la G y la A). En el laboratorio se usan como «tijeras» para cortar y pegar fragmentos de ADN. Según dónde corten, dejan **extremos cohesivos** (con un trozo de hebra sencilla que sobresale, que facilita pegar fragmentos) o **extremos romos**.

Los **códigos IUPAC** permiten buscar motivos con posiciones variables: R = A/G, Y = C/T, W = A/T, S = C/G, K = G/T, M = A/C y N = cualquier base. Por ejemplo, la caja TATA se busca como `TATAWAW`.

#### Cómo usarla, paso a paso

1. Selecciona la secuencia y abre **Motif search**.
2. Escribe en **Motif or site name** una secuencia (`GAATTC`), un motivo con códigos IUPAC (`TATAWAW`) o el nombre de una enzima (`EcoRI`).
3. Pulsa **Search** (Buscar), o haz clic en uno de los botones rápidos que aparecen debajo con las enzimas y motivos más habituales.

#### Cómo leer los resultados

- Un resumen con el número de coincidencias en la hebra + y en la −.
- Un **mapa** con una marca en cada coincidencia.
- Una **tabla** con la posición, la hebra y la secuencia exacta encontrada (útil con códigos IUPAC). Si buscaste una enzima, una columna **Cut** indica entre qué dos bases corta, y una nota explica el tipo de extremo que deja.
- Enzimas incluidas: EcoRI, BamHI, HindIII, NotI, EcoRV, XhoI, PstI, SmaI, KpnI, SacI, XbaI, SalI, NcoI, NdeI SpeI y BglII; además del motivo TATA-box.
- La mayoría de los sitios de restricción son **palíndromos** (se leen igual en las dos hebras, como GAATTC), así que cada sitio aparece una vez en la hebra + y otra en la −: el número de coincidencias es el doble del número de sitios.

> **Consejo:** Para clonar un gen en un plásmido conviene elegir enzimas que **no** corten dentro del gen. Busca cada candidata y descarta las que tengan coincidencias en él.

### Pestaña «K-mers» (Palabras de longitud k)

#### ¿Qué hace?

Cuenta todas las «palabras» de **k** letras seguidas que aparecen en la secuencia (solapándose) y muestra las más frecuentes. Por ejemplo, con k = 3 en `ATGCA` las palabras son ATG, TGC y GCA.

#### ¿Para qué sirve?

- Detectar **repeticiones** y secuencias que aparecen más de lo esperado.
- Comparar secuencias o identificar organismos por su «firma» de k-mers (muchas herramientas modernas de ensamblaje y clasificación de genomas se basan en esto).
- Ver la diversidad de una secuencia: cuántas palabras distintas tiene frente a las posibles.

#### Cómo usarla, paso a paso

1. Selecciona la secuencia y abre **K-mers**.
2. Elige **k-mer length** (longitud k, de 1 a 20). Con k = 1 obtienes la frecuencia de cada base; con k = 2, de cada pareja (dinucleótidos), etc.
3. Elige en **Show top** cuántas palabras quieres ver (las más frecuentes).

#### Cómo leer los resultados

- El total de k-mers, cuántos son distintos y cuántos serían posibles (4 elevado a k).
- Una tabla ordenada por frecuencia con el k-mer, su número de apariciones y su porcentaje, con una barra.
- Descarga en **CSV**.

### Pestaña «Composition» (Composición)

#### ¿Qué hace?

Analiza cómo cambia la composición de la secuencia **a lo largo de ella**, dividiéndola en ventanas que se desplazan, y busca **islas CpG**.

#### ¿Para qué sirve?

- Ver regiones ricas o pobres en GC (pueden corresponder a genes, a ADN de otra procedencia o a zonas difíciles de secuenciar o amplificar por PCR).
- Localizar el **origen de replicación** en genomas de bacterias y virus mediante el sesgo GC.
- Detectar **repeticiones y zonas de baja complejidad**.
- Encontrar **islas CpG**, que en vertebrados suelen marcar el inicio de los genes.

#### Un poco de biología

**Ventana deslizante:** en lugar de un solo valor para toda la secuencia, se calcula el valor en un trozo (la ventana, por ejemplo de 100 bp), luego se desplaza un poco y se vuelve a calcular. Así se obtiene un gráfico a lo largo de la secuencia.

**Sesgo GC (GC skew)** = (G − C) / (G + C). Mide si en una hebra hay más G que C. Varía entre −1 y 1. En muchas bacterias y fagos cambia de signo cerca del origen y del final de la replicación, porque las dos hebras se copian de forma distinta. El **sesgo AT** es lo mismo con A y T.

**Complejidad (entropía de Shannon):** mide cuán variada es la composición, en bits. Vale 2 cuando las cuatro bases aparecen por igual y baja hacia 0 en zonas repetitivas (por ejemplo, AAAAAAA o ATATATAT).

**Islas CpG:** «CpG» es una C seguida de una G en la misma hebra. En el ADN de vertebrados este dinucleótido es raro (su citosina suele estar metilada y muta con facilidad), salvo en zonas llamadas islas CpG, que con frecuencia están en los **promotores** de los genes. La web usa los criterios clásicos de Gardiner-Garden y Frommer: al menos **200 bp**, **≥ 50 % de GC** y una relación **CpG observado/esperado ≥ 0,6**. En bacterias y fagos el CpG no escasea, así que salen muchas «islas» que no tienen significado biológico.

#### Cómo usarla, paso a paso

1. Selecciona la secuencia y abre **Composition**.
2. Ajusta **Window size** (Tamaño de ventana) en bp. Ventanas pequeñas muestran más detalle; ventanas grandes suavizan el ruido en secuencias largas.

#### Cómo leer los resultados

- Recuadros con los valores globales: sesgo GC, sesgo AT, complejidad (máximo 2 bits), relación CpG observado/esperado (en ADN humano típico ronda 0,2–0,3) y número de islas CpG.
- Tres gráficos a lo largo de la secuencia: **% GC** (con las islas CpG sombreadas), **sesgo GC** y **complejidad** (las caídas marcan repeticiones).
- Si hay islas CpG, una tabla con su posición, longitud, % GC y relación observado/esperado.
- Descarga de los valores de cada ventana en **CSV**.

### Pestaña «Tm & weight» (Temperatura de fusión y peso)

#### ¿Qué hace?

Calcula la **temperatura de fusión (Tm)** de la secuencia y su **peso molecular**, como ADN de una hebra, de doble hebra y el de la proteína del ORF más largo.

#### ¿Para qué sirve?

- **Diseñar cebadores (primers) para PCR:** los dos cebadores de una pareja deben tener una Tm parecida, normalmente entre 52 y 65 °C, y la temperatura de hibridación de la PCR se fija unos grados por debajo.
- Convertir entre masa y número de moléculas (por ejemplo, para preparar diluciones de ADN o calcular moles de un oligonucleótido).
- Estimar el tamaño de la proteína esperada (por ejemplo, para saber en qué banda aparecerá en un gel de electroforesis).

#### Un poco de biología

La **Tm** es la temperatura a la que la mitad de las moléculas de ADN de doble hebra se han separado en hebras sencillas. Depende de la secuencia (los pares G–C aportan más estabilidad), de la concentración de sal (el Na⁺ estabiliza la doble hélice) y de la concentración del cebador.

**Método del vecino más próximo (nearest-neighbor):** es el que usan los fabricantes de cebadores. Suma la energía de cada pareja de bases vecinas (parámetros de SantaLucia, 1998) y corrige según la sal y la concentración de cebador. Es el valor recomendado.

**Fórmula básica:** una estimación rápida que ignora la sal. Por debajo de 14 nt usa la regla de Wallace, 2·(A+T) + 4·(G+C); desde 14 nt, 64,9 + 41·(G+C − 16,4)/longitud.

Ambos cálculos están pensados para **oligonucleótidos de hasta unos 60 nt**. Para secuencias más largas la web avisa de que el valor no es significativo.

#### Cómo usarla, paso a paso

1. Carga la secuencia del cebador (por ejemplo, pégala sola en el cuadro) y abre **Tm & weight**.
2. Ajusta **Na⁺** (concentración de sodio, en mM) y **Primer** (concentración del cebador, en nM) según tu protocolo. Los valores por defecto (50 mM y 250 nM) son los típicos de una PCR.
3. El campo de **longitud mínima del ORF** solo afecta al recuadro del peso de la proteína.

#### Cómo leer los resultados

- **Tm, nearest-neighbor:** la Tm recomendada, en °C (no se calcula si la secuencia tiene N o menos de 2 bases; la web explica el motivo).
- **Tm, basic formula:** la estimación rápida, indicando qué fórmula se usó.
- **DNA weight, single / double strand:** peso del ADN de una hebra y de doble hebra, en kDa y g/mol.
- **Longest ORF protein:** peso de la proteína del ORF más largo, en kDa, con su longitud y posición.

> **Consejo:** El archivo `data/samples/primers.fasta` contiene cebadores de ejemplo para practicar con esta pestaña.

### Pestaña «Read quality» (Calidad de las lecturas)

#### ¿Qué hace?

Analiza la **calidad** de las lecturas de un archivo **FASTQ** procedente de un secuenciador, y permite **recortar** las partes de baja calidad y descargar las lecturas limpias. Esta pestaña solo funciona con archivos FASTQ; con un FASTA muestra un aviso.

Es la pestaña que responde a la pregunta «¿son fiables mis datos de secuenciación?».

#### ¿Para qué sirve?

- Comprobar si una secuenciación ha salido bien **antes** de usar los datos en otros análisis (ensamblaje, búsqueda de mutaciones, etc.). Datos de mala calidad producen resultados erróneos.
- Ver **dónde** baja la calidad dentro de las lecturas (normalmente al final).
- **Limpiar** los datos: cortar los extremos malos y descartar las lecturas que quedan demasiado cortas.

#### Un poco de biología

**¿Qué es la secuenciación?** Las máquinas de secuenciación (como las de Illumina) no leen un genoma de una sola vez: lo trocean en millones de fragmentos y leen cada fragmento, base a base, para producir **lecturas** («reads») de, por ejemplo, 100–300 bases. En cada posición la máquina mide una señal (normalmente de luz) y decide qué base es. A veces la señal es dudosa, y la máquina puede equivocarse.

**¿Cómo se mide la calidad?** Por eso, junto a cada base, la máquina anota cuán segura está mediante la **puntuación de calidad Phred (Q)**. Se define como Q = −10 · log₁₀(P), donde P es la probabilidad de que la base esté mal leída. Dicho de forma práctica: **una puntuación Q significa una probabilidad de error de 1 entre 10 elevado a Q/10**. Cada 10 puntos más, la probabilidad de error es 10 veces menor:

| Calidad Phred | Probabilidad de error | Precisión | Interpretación |
|---|---|---|---|
| Q10 | 1 de cada 10 | 90 % | Mala |
| Q20 | 1 de cada 100 | 99 % | Aceptable (umbral mínimo habitual) |
| Q30 | 1 de cada 1 000 | 99,9 % | Buena (estándar de referencia en Illumina) |
| Q40 | 1 de cada 10 000 | 99,99 % | Excelente |

En el archivo FASTQ, cada puntuación se guarda como un carácter (codificación Phred+33: el carácter `I` equivale a Q40, `5` a Q20, `+` a Q10). La web hace esta conversión por ti.

**¿Por qué baja la calidad al final de las lecturas?** En la secuenciación por síntesis, en cada ciclo algunas moléculas de un mismo grupo se adelantan o se retrasan y la señal se va mezclando y debilitando. Por eso es normal que las lecturas de Illumina pierdan calidad hacia su extremo 3′ (el final).

**¿Qué es el recorte (trimming)?** Consiste en cortar la parte final de baja calidad de cada lectura. La web usa el método de **ventana deslizante** del programa Trimmomatic (SLIDINGWINDOW): recorre la lectura desde el principio con una ventana de 4 bases y, en cuanto la calidad media de la ventana cae por debajo del umbral (Q20), corta la lectura ahí (conservando las bases iniciales de esa ventana que sí superan el umbral). Después descarta las lecturas que han quedado más cortas que el mínimo (36 bp), porque son demasiado cortas para ser útiles.

#### Cómo usarla, paso a paso

1. Carga un archivo **FASTQ** (botón **Upload file** o arrastrándolo). La web abre automáticamente esta pestaña. Para probar, usa `data/samples/lambda_reads.fastq`.
2. Revisa los recuadros del resumen y los dos gráficos.
3. En los controles de recorte ajusta: **Trim below** (recortar por debajo de Q…, por defecto Q20), **Window** (ventana, por defecto 4 bases) y **Keep reads of at least** (conservar lecturas de al menos… bp, por defecto 36).
4. Revisa los recuadros «después del recorte» para ver cuántos datos se conservan y cuánto mejora la calidad.
5. Pulsa **Trimmed reads (FASTQ)** para descargar las lecturas recortadas, listas para otros programas.

#### Cómo leer los resultados

- **Reads:** número de lecturas cargadas.
- **Read length:** longitud de las lecturas (o el rango, si varían, con su media).
- **Mean quality:** calidad media de todas las bases. Por encima de Q30 es muy buena.
- **≥ Q30** y **≥ Q20:** porcentaje de bases con calidad igual o superior a 30 (error de 1 en 1 000) y a 20 (error de 1 en 100). Es el indicador más usado: una buena secuenciación Illumina suele tener más del 75–80 % de bases ≥ Q30.
- **GC · N:** porcentaje de GC de todas las lecturas y porcentaje de bases N (no leídas). Un % de N alto indica problemas.
- **Gráfico «Quality along the reads» (Calidad a lo largo de las lecturas):** al estilo del programa FastQC. Para cada posición dentro de la lectura muestra una caja con la distribución de calidades (la caja va del primer al tercer cuartil, la línea central es la mediana, los bigotes van del percentil 10 al 90). El fondo está coloreado: **verde** = buena (Q28 o más), **ámbar** = aceptable, **rojo** = mala (por debajo de Q20). Si las cajas bajan a la zona roja hacia la derecha, conviene recortar. En lecturas largas, las posiciones se agrupan en bloques.
- **Histograma «Reads by mean quality» (Lecturas por calidad media):** cuántas lecturas tienen cada calidad media. Lo ideal es un solo pico en la zona verde; una cola en rojo indica un grupo de lecturas malas.
- **Quality trimming (Recorte por calidad):** recuadros con las lecturas conservadas, el % de bases conservadas, la longitud media tras el recorte, y el % ≥ Q30 y la calidad media **después**, comparados con los valores de antes.
- Descargas: **Per-position CSV** (los datos del gráfico de cajas) y **Trimmed reads (FASTQ)** (las lecturas limpias).

> **Consejo:** Si al recortar se pierde un porcentaje muy alto de lecturas, el problema no es el recorte sino la secuenciación. Si se pierde poco y el % ≥ Q30 sube, el recorte ha cumplido su función.

> **Además:** Cada lectura del FASTQ se carga también como secuencia en la lista de la izquierda, así que puedes seleccionar una lectura concreta y analizarla con las demás pestañas.

### Pestaña «All sequences» (Todas las secuencias)

#### ¿Qué hace?

Muestra una tabla resumen con **todas** las secuencias cargadas a la vez (no solo la seleccionada).

#### ¿Para qué sirve?

- Comparar varias secuencias de un vistazo (longitud, GC, número de genes).
- Obtener un informe completo para guardar o compartir. Es el mismo informe que produce la versión de terminal con `python3 main.py archivo.fasta`.

#### Cómo usarla, paso a paso

1. Carga un archivo con varias secuencias (por ejemplo, `data/samples/human_genes.fasta`, con los ARNm humanos de la insulina y de p53).
2. Abre **All sequences**.
3. Ajusta la **longitud mínima del ORF** si quieres cambiar qué ORF se cuentan.

#### Cómo leer los resultados

- Una fila por secuencia con: nombre, longitud, % GC, número de A, C, G y T, número de ORF encontrados y el **ORF más largo** (posición, hebra y longitud de la proteína).
- Descarga de toda la tabla en **CSV**.

## 5. Ejemplos de uso de principio a fin

### 5.1 Analizar un gen desconocido

1. Pega la secuencia y pulsa **Analyze**.
2. En **Overview** comprueba la longitud y el % de GC.
3. En **ORFs** busca el ORF más largo: probablemente es el gen. Anota su posición y su hebra.
4. Descarga la proteína en FASTA y búscala en BLAST para saber qué es.
5. En **Translation** comprueba la proteína en el marco correspondiente.

### 5.2 Preparar una clonación

1. Carga la secuencia del gen.
2. En **Motif search** prueba las enzimas candidatas y descarta las que cortan dentro del gen.
3. Diseña los cebadores con las enzimas elegidas y comprueba su temperatura en **Tm & weight**: los dos deben tener una Tm parecida.

### 5.3 Revisar los datos de una secuenciación

1. Sube el archivo FASTQ (por ejemplo, `data/samples/lambda_reads.fastq`).
2. En **Read quality** mira el % ≥ Q30 y el gráfico de calidad por posición.
3. Si la calidad baja al final, ajusta el recorte y descarga las lecturas recortadas.
4. Usa ese archivo limpio en los siguientes pasos de tu análisis.

### 5.4 Archivos de ejemplo incluidos

| Archivo (en data/samples/) | Úsalo para probar… |
|---|---|
| `lambda_phage.fasta` | Genoma del fago lambda: ORFs, composición (sesgo GC), motivos. |
| `lambda_reads.fastq` | Lecturas de secuenciación: pestaña Read quality y recorte. |
| `human_genes.fasta` | ARNm humanos de la insulina y de p53: All sequences, ORFs, peso de proteínas. |
| `primers.fasta` | Cebadores cortos: Tm & weight. |
| `rna_example.fasta` | Secuencias de ARN (con U). |
| `feature_tour.fasta` | Un recorrido por varias funciones de la web. |
| `genbank_style_paste.txt` | Secuencia con números y espacios, como se copia de GenBank. |

## 6. Glosario

| Término | Significado |
|---|---|
| Aminoácido | Cada una de las 20 piezas que forman las proteínas. |
| bp / nt | Pares de bases / nucleótidos: unidades de longitud del ADN. |
| Cebador (primer) | Fragmento corto de ADN (unos 18–30 nt) que marca dónde empieza la copia en una PCR. |
| Codón | Grupo de tres bases que corresponde a un aminoácido o a una parada. |
| Complementaria inversa | La otra hebra del ADN, leída de 5′ a 3′. |
| Contenido en GC | Porcentaje de bases G y C de una secuencia. |
| CpG / isla CpG | Una C seguida de una G; región rica en CpG, a menudo en los promotores de los genes de vertebrados. |
| Enzima de restricción | Proteína que corta el ADN en una secuencia concreta. |
| FASTA / FASTQ | Formatos de archivo de secuencias; FASTQ incluye la calidad de cada base. |
| k-mer | Subsecuencia de k bases seguidas. |
| Lectura (read) | Fragmento de ADN leído por un secuenciador. |
| Marco de lectura | Cada una de las tres formas (seis con la otra hebra) de dividir una secuencia en codones. |
| Motivo | Secuencia corta con significado biológico que se busca dentro de otra. |
| ORF | Marco abierto de lectura: de un ATG a un codón de parada; un posible gen. |
| PCR | Técnica para hacer millones de copias de un fragmento de ADN. |
| Phred (Q) | Puntuación de calidad de una base: Q = −10·log₁₀(probabilidad de error). |
| Tm | Temperatura de fusión: aquella a la que la mitad del ADN de doble hebra se ha separado. |
| Transcripción / traducción | Paso de ADN a ARN / de ARN a proteína. |
| Trimming (recorte) | Eliminación de las partes de baja calidad de las lecturas. |
