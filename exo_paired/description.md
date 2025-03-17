# Обработка экзома  
  

---

## Задача
**Дано**:
- Парный экзом: норма/опухоль, использумеый экзомный набор: NEXome Plus Panel v1.0
- Образец: пациентка с тройным негативным раком молочной железы (ТНРМЖ/TNBC)

**Цель**: определение герминальных и соматических вариантов, предположение терапии

---

##  Ход работы  
1. QC по образцам: образцы высокого качества с особенностями, характерными для экзомов. Также есть достаточно большое количество дупликатов.

2. Получены варианты с помощью [скрипта](norm_vcf.sh) и проаннотированы VEP+wANNOVAR, отфильтрованы [R-скриптом](anotation.R)
```
chr	start	end	ref	alt	gene_ref_gene_with_ver	exonic_func_ref_gene_with_ver	otherinfo6
chr3	119582455	119582455	G	-	ADPRH	frameshift deletion	NA
chr5	173609434	173609434	-	TGAACCACACTCTGCCTCAAACCATTTCGCAATTGGTTTTTGTTCATTGTAGGATTCCATTCCTG	BOD1	stopgain	NA
chr6	169602251	169602251	G	A	WDR27	stopgain	rs151208349
chr12	57525783	57525784	TG	-	MBD6	frameshift deletion	NA
chr19	4511759	4511759	-	CA	PLIN4	frameshift insertion	NA
chr19	4511762	4511763	AC	-	PLIN4	frameshift deletion	NA
chr19	57641161	57641161	A	-	ZNF211	frameshift deletion	NA
chr20	50594589	50594590	TG	-	RIPOR3	frameshift deletion	NA
chr22	26496858	26496858	G	A	TFIP11	stopgain	rs374892403
```

Из найденных вариантов наибольший интерес представляет вставка в третьем экзоне гена BOD1, ответственном за ориентацию сестринских хроматид в митозе, что может быть связано с образованием опухоли путём формирования хромосомной нестабильности. Вариант не встречается в БД и при использовании Varsome/Franklin. Вставка была обработана через blast со следующим результатом:
```
>XM_047443198.1:28-92 PREDICTED: Homo sapiens biorientation of chromosomes in cell division protein 1 (LOC124905450), mRNA
CAGGAATGGAATCCTACAATGAACAAAAACCAATTGCGAAATGGTTTGAGGCAGAGTGTGGTTCA
>AC208950.3:16625-16689 Homo sapiens FOSMID clone ABC10-43645600F20 from chromosome 18, complete sequence
CAGGAATGGAATCCTACAATGAACAAAAACCAATTGCGAAATGGTTTGAGGCAGAGTGTGGTTCA
```


3. Колинг соматических вариантов с помощью Mutect2 и фильтрация проводились по [скрипту](tumor_to_vcf.sh). Подбор параметров для фильтрации проходил с использование [небольшого скрипта](mean.sh)

4. Соматические варианты проаннотированы через VEP, отфильтрованы тем же R-скриптом.
```
location	allele	consequence	symbol	exon	c_dna_position	existing_variation	ref_allele	uploaded_allele
12:121240524-121240526	-	frameshift_variant	CAMKK2	16/16	1618-1619	rs63023660	TT	TT/-/T/TTT
12:121240524-121240526	-	frameshift_variant	CAMKK2	11/11	1041-1042	rs63023660	TT	TT/-/T/TTT
2:74958681-74958684	-	start_lost,inframe_deletion,start_retained_variant	POLE4	1/4	40-42	rs528770146	GGC	GGC/-
4:5989141-5989145	-	frameshift_variant	C4orf50	28/34	8297-8300	rs143859826	CTCT	CTCT/-/CT/CTCTCT
4:5989141-5989145	CTCTCT	frameshift_variant	C4orf50	28/34	8297-8300	rs143859826	CTCT	CTCT/-/CT/CTCTCT
4:5989141-5989145	-	frameshift_variant,NMD_transcript_variant	C4orf50	1/8	920-923	rs143859826	CTCT	CTCT/-/CT/CTCTCT
4:5989141-5989145	CTCTCT	frameshift_variant,NMD_transcript_variant	C4orf50	1/8	920-923	rs143859826	CTCT	CTCT/-/CT/CTCTCT
4:5989141-5989145	-	frameshift_variant	C4orf50	6/12	2901-2904	rs143859826	CTCT	CTCT/-/CT/CTCTCT
4:5989141-5989145	CTCTCT	frameshift_variant	C4orf50	6/12	2901-2904	rs143859826	CTCT	CTCT/-/CT/CTCTCT
```
Полученные варианты не удалось однозначно проинтерпретировать. Возможные причины:
- слишком строгая фильтрация
- использование неподходящих средств аннотации

# Результат
\**в работе*\* 