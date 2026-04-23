# Plan testów E2E dla annotacji pluginami

Ten dokument opisuje plan walidacji end-to-end dla pięciu cache pluginów wspieranych obecnie przez lokalny workflow `vepyr`:

- `clinvar`
- `spliceai`
- `cadd`
- `alphamissense`
- `dbnsfp`

Celem jest potwierdzenie, że każdy plugin można zbudować z realnego lokalnego pliku źródłowego, zmaterializować w oczekiwanym layoucie cache, użyć w annotacji `vepyr`, a następnie porównać wynik z Ensembl VEP przez `vepyr-diffly`.

Pełne generowanie cache może być kosztowne, szczególnie dla CADD, SpliceAI i dbNSFP. W trakcie rozwoju używamy uruchomień ograniczonych do chromosomu oraz `--preview-rows`, żeby walidować poprawność i mierzyć względną wydajność. Finalna akceptacja powinna bazować na pełnych wejściach pluginów tam, gdzie jest to wykonalne, albo na jasno opisanym ograniczonym zakresie, jeśli pełny run jest zbyt duży dla lokalnej maszyny.

## Procedura testowa

1. Użyć istniejącego pełnego core cache parquet jako stałej bazy annotacji.
2. Zbudować pełny cache dla jednego pluginu naraz z lokalnego pliku źródłowego, na przykład tylko `clinvar`.
3. Zapisać rozmiar skompresowanego wejścia, rozmiar wygenerowanego parquet, rozmiar wygenerowanego Fjall oraz czas budowy cache.
4. Uruchomić annotację Ensembl VEP na tych samych wariantach z włączonym odpowiadającym pluginem.
5. Uruchomić annotację `vepyr` z wygenerowanym lokalnym cache i tym samym pluginem.
6. Porównać oba annotowane VCF-y przez `vepyr-diffly`, koncentrując akceptację na zgodności danych pluginowych.
7. Oznaczyć plugin jako zgodny z VEP tylko wtedy, gdy semantyczne porównanie raportuje 100% zgodności danych pluginowych dla stosunkowo dużej próbki.
8. Dopiero po osiągnięciu 100% zgodności zapisać finalne dane pomiarowe w tabeli wyników i oznaczyć plugin jako działający.

Porównanie musi zapisywać oba całkowite czasy annotacji:

- Czas VEP to pełny wall-clock time kroku annotacji Ensembl VEP.
- Czas `vepyr` to pełny wall-clock time lokalnej annotacji `vepyr` z użyciem wygenerowanego cache.

## Fixture VCF

Pierwszy większy fixture dla testu ClinVar:

- Ścieżka: `/tmp/clinvar_e2e_stratified_1m.vcf`
- Źródło: `plugins/clinvar.vcf.gz`
- Strategia: deterministyczna próbka stratifikowana, maksymalnie `40_000` rekordów na chromosom/kontig.
- Rozmiar i liczność: `923_845` rekordów, `45` linii nagłówka, `28M`.
- Przygotowanie: kolumna `INFO` została wyczyszczona do `.`, żeby wejściowy ClinVar nie był mylony z annotacją dodaną przez VEP/`vepyr`.
- Walidacja: podstawowy check VCF potwierdził `0` rekordów z błędną liczbą kolumn albo niewyczyszczonym `INFO`.

Golden sample z GIAB HG002 do realistycznych benchmarków:

- Ścieżka: `/tmp/hg002_grch38_100k_stratified.vcf`
- Źródło: `/Users/lukaszjezapkowicz/Downloads/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf`
- Strategia: deterministyczna próbka stratifikowana po chromosomach 1-22; `4_546` rekordów dla chromosomów 1-10 i `4_545` dla 11-22.
- Rozmiar i liczność: `100_000` rekordów, `230` linii nagłówka, `67M`.
- Przygotowanie: pełne kolumny `INFO`, `FORMAT` i próbka `HG002` są zachowane; nazwy chromosomów znormalizowano z `chrN` do `N`, żeby plik był gotowy pod lokalne VEP/cache GRCh38.
- Walidacja: `0` rekordów z błędną liczbą kolumn i `0` pozostałych rekordów z prefiksem `chr`; skład alleli: `87_273` SNV, `6_672` insercji, `6_880` delecji, `9` MNP.

## Tabela wyników

| Plugin | Zakres testu | Czas tworzenia cache | Rozmiar wejścia pluginu | Rozmiar outputu parquet | Rozmiar outputu Fjall | 100% zgodność z VEP | Czas annotacji VEP | Czas annotacji vepyr | Uwagi |
| --- | --- | ---: | ---: | ---: | ---: | :---: | ---: | ---: | --- |
| `clinvar` | Full cache + plugin E2E on `/tmp/clinvar_e2e_stratified_100k.vcf` | 640.7s convert / 648.4s wall | 175.4MB gz / 1.8GB raw | 165M | 760M | tak (plugin-only compare, 100k) | 2470.36s real / 2395.16s user / 45.80s sys | 397.39s real / 306.91s user / 79.41s sys | 2026-04-19: plugin-only compare na `runs/clinvar-e2e-100k-bench/compare-plugins-only` dał `variant.equal=true`, `consequence.equal=true`, `99988` joined-equal na obu tierach, `0/0/0` mismatchy. Outputy annotacji: oba VCF-y ~`1.3G`, po `100042` linii. Full compare `--everything` nadal pokazuje bazowy drift poza pluginem (`compare/summary.json`), więc akceptacja dotyczy danych pluginowych zgodnie z ustaloną procedurą. Wrapper `/usr/bin/time -l` zwrócił błąd sandboxowy przy `sysctl kern.clockrate`, więc peak RSS nie został zapisany. |
| `spliceai` | Full cache parquet-only + plugin E2E on `/tmp/hg002_grch38_100k_stratified.vcf` | 7,484.2s convert / 7,484.2s wall | 26.6GB gz / raw `UNKNOWN` (`gzip -l` unreliable for BGZF/multi-member gzip) | 87G | n/a (`--no-plugin-fjall`) | tak (plugin-only compare, 100k) | ~1532s wall (`2026-04-22T10:52:29+02:00` -> `2026-04-22T11:18:01+02:00`) | ~849s wall (`2026-04-22T11:41:26+02:00` -> `2026-04-22T11:55:35+02:00`, vepyr-only rerun after symbol-gating fix) | 2026-04-22: full cache has `3,393,685,728` rows in `24` parquet files (`93,687,625,027` bytes by parquet file metadata). Initial full-cache runtime smoke exposed that parquet plugin lookup was loading the whole per-chromosome SpliceAI file (`chr1.parquet=7.7G`) for one input variant; runtime now collects target `(pos, ref, alt)` keys from the input VCF and filters plugin parquet reads to those keys. The first 100k run had plugin-only consequence drift (`203` left-only, `741` right-only) because SpliceAI payload was attached by variant only; VEP gates SpliceAI rows by consequence gene symbol. `vepyr` now selects SpliceAI candidate rows by the current consequence `SYMBOL` and emits empty plugin fields for non-matching transcript/regulatory rows. Fixed compare `runs/spliceai-e2e-100k-bench/compare-plugins-only-symbolfix/summary.json` passed fully: variant tier equal (`100834` joined-equal), plugin consequence tier equal (`36224` joined-equal), `0/0/0` mismatches. Outputy annotacji: `vep.annotated.vcf=606M` (`101076` linii), fixed `vepyr.annotated.symbolfix.vcf=599M` (`101063` linii). Compare timings: `variant_summary=12.579s`, `consequence_bucketization=11.214s`, `variant_diff=0.458s`, `consequence_diff=0.786s`. Positive 10-variant smoke on `/tmp/spliceai_positive_10.vcf` also passed. Peak RSS: `UNCONFIRMED`. |
| `cadd` | chr1-only parquet cache + plugin E2E on `/tmp/hg002_chr1_cadd_100k.vcf` | chr1 completed during interrupted full build; exact wall `UNCONFIRMED` | 81G SNV + 1.2G indel (`.tbi`: 2.6M + 1.8M) | 21G (`cadd/chr1.parquet`, `699,768,898` rows) | n/a (`--no-plugin-fjall`) | tak (plugin-only compare, chr1 100k) | ~1272s wall (`2026-04-23T12:55:25+02:00` -> `2026-04-23T13:16:37+02:00`) | ~156s wall (`2026-04-23T13:46:37+02:00` -> `2026-04-23T13:49:13+02:00`, vepyr-only rerun after CADD SNV-fallback fix) | 2026-04-23: CADD full cache generation was stopped after the complete `chr1.parquet` write; `chr10.parquet` may be partial and is not part of this result. Upstream CADD `.tbi` files were downloaded instead of locally indexing the 81G SNV source. Initial chr1 100k compare had `19` VEP-only plugin consequence rows, all equal-length alleles with one changed base such as `GATT>TATT`; VEP's CADD plugin resolves these as single SNV lookups, while vepyr originally required exact full-allele lookup. Runtime now adds a single-base-substitution fallback and includes those fallback keys in parquet prefiltering, without rebuilding cache. Fixed compare `runs/cadd-chr1-e2e-100k/compare-plugins-only-caddfix/summary.json` passed fully: variant tier equal (`101350` joined-equal), plugin consequence tier equal (`100450` joined-equal), `0/0/0` mismatches. Outputy annotacji: `vep.annotated.vcf=483M` (`101585` linii), fixed `vepyr.annotated.caddfix.vcf=482M` (`101579` linii). Compare timings: `variant_summary=9.995s`, `consequence_bucketization=8.739s`, `variant_diff=0.519s`, `consequence_diff=0.761s`. Smoke 100 on `/tmp/hg002_chr1_cadd_100.vcf` also passed (`100` variant joined-equal, `99` plugin consequence joined-equal). Peak RSS: `UNCONFIRMED`; `/usr/bin/time -l` is blocked by `sysctl kern.clockrate` in this environment. |
| `alphamissense` | Full cache parquet-only + plugin E2E on `/tmp/hg002_grch38_100k_stratified.vcf` | 661.0s convert / 661.0s wall | 591.4MB gz / 1.1GB raw | 2.0G | n/a (`--no-plugin-fjall`) | tak (plugin-only compare, 100k) | ~896s wall (`2026-04-20T08:26:17+02:00` -> `2026-04-20T08:41:13+02:00`) | ~759s wall (`2026-04-20T08:41:13+02:00` -> `2026-04-20T08:53:52+02:00`) | 2026-04-20: `runs/alphamissense-e2e-100k-bench/summary.json` dał `variant.equal=true`, `consequence.equal=true`, `427` joined-equal consequence rows, `0/0/0` mismatchy przy `--compare-only-plugins`. Outputy annotacji: `vep.annotated.vcf=579M` (`101069` linii), `vepyr.annotated.vcf=585M` (`101063` linii). Dodatni smoke fixture `/tmp/alphamissense_positive_10.clean.vcf` też dał `10/10` joined-equal. Peak RSS: `UNCONFIRMED`. |
| `dbnsfp` | Full cache parquet-only + plugin E2E on `/tmp/hg002_grch38_100k_stratified.vcf` | 36,754.8s convert / 36,754.8s wall | 46.8GB gz / raw `UNKNOWN` (`gzip -l` unreliable for BGZF/multi-member gzip) | 5.7G | n/a (`--no-plugin-fjall`) | tak (plugin-only compare, 100k) | ~940s wall (`2026-04-21T20:20:16+02:00` -> `2026-04-21T20:35:56+02:00`) | ~928s wall (`2026-04-21T20:35:56+02:00` -> `2026-04-21T20:51:24+02:00`) | 2026-04-21: fixed full 100k run `runs/dbnsfp-e2e-100k-bench-fixed` passed under `--compare-only-plugins`: variant tier equal (`100834` joined-equal), plugin consequence tier equal (`613` joined-equal), `0/0/0` mismatches; `consequence_mismatches.tsv` contains only the header. Outputy annotacji: `vep.annotated.vcf=606M` (`101086` linii), `vepyr.annotated.vcf=615M` (`101063` linii). Compare timings: `variant_summary=12.814s`, `consequence_bucketization=11.585s`, `variant_diff=0.538s`, `consequence_diff=0.118s`. Cache has `86,435,706` rows in `25` parquet files. Positive smoke `/tmp/dbnsfp_positive_100.vcf` also passed (`variant=yes`, `consequence=yes`, `88` joined-equal, `0/0/0`), and mini regression `runs/dbnsfp-revel-probe-mini-fixed` passed after fixing VEP-style string preservation, dbNSFP separator escaping, and per-consequence duplicate-row selection. Peak RSS: exact value `UNCONFIRMED`; sampled `%MEM` reached about VEP `6.3%` and `vepyr` `29.4%`, so dbNSFP vepyr runtime is notably memory-heavier than smaller plugins. |

## Baseline bez pluginów

Ten baseline służy tylko jako punkt odniesienia dla testów pluginowych na tym samym golden sample. Nie jest liczony jako wynik żadnego pluginu.

| Zakres testu | Core cache | Zgodność VEP vs vepyr | Czas annotacji VEP | Czas annotacji vepyr | Czas compare | Rozmiar outputu VEP | Rozmiar outputu vepyr | Pamięć | Uwagi |
| --- | ---: | :---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| Root cache only on `/tmp/hg002_grch38_100k_stratified.vcf` | 43G | variant: tak; consequence: nie (`19/19/0`, `1_422_639` joined-equal) | ~876s wall (`2026-04-20T12:42:27+02:00` -> `2026-04-20T12:57:03+02:00`) | ~700s wall (`2026-04-20T12:57:03+02:00` -> `2026-04-20T13:08:43+02:00`) | ~154s wall; summary timings: `11.585s + 34.920s + 0.036s + 0.332s + 106.861s` | 577M (`101067` linii) | 584M (`101063` linii) | Peak RSS: `UNCONFIRMED`; sampled `%MEM`: VEP ~8.3%, vepyr ~9.6% | 2026-04-20: `runs/root-only-e2e-100k-bench-rerun/summary.json`. Run bez pluginów na tym samym 100k HG002 goldenie. Variant tier ma 100% zgodności (`100834` joined-equal). Mały baseline drift consequence jest poza oceną plugin-only E2E, ale warto go uwzględniać przy interpretacji full compare. |

## Zasady pomiaru

- `Zakres testu` musi mówić, czy wynik dotyczy pełnego wejścia, konkretnego chromosomu, `--preview-rows`, albo kombinacji, na przykład `chr1 preview_rows=1_000_000`.
- `Czas tworzenia cache` powinien obejmować przygotowanie/wycinanie źródła oraz generowanie parquet i Fjall.
- `Rozmiar wejścia pluginu` powinien oznaczać rozmiar lokalnego skompresowanego źródła, a dla pluginów wieloplikowych takich jak CADD sumę skompresowanych źródeł.
- `Rozmiar outputu parquet` powinien oznaczać łączny rozmiar wygenerowanego katalogu parquet pluginu.
- `Rozmiar outputu Fjall` powinien oznaczać łączny rozmiar wygenerowanego katalogu `<plugin>.fjall`.
- `100% zgodność z VEP` ustawiamy na `tak` tylko wtedy, gdy `vepyr-diffly` raportuje pełną zgodność dla wybranego zakresu; w przeciwnym razie wpisujemy `nie` i linkujemy artefakt z mismatchami w `Uwagi`.
- Dla akceptacji pluginu najważniejsza jest zgodność pól pluginowych; różnice w niepowiązanych polach bazowej annotacji należy opisać osobno i nie mieszać z wynikiem zgodności pluginu.
- Jeśli pełny run nie jest lokalnie wykonalny, zostawiamy wynik dla ograniczonego zakresu, ale zakres musi być jawnie opisany.
- Jeżeli środowisko blokuje peak RSS z `/usr/bin/time -l` (`sysctl kern.clockrate: Operation not permitted`), zapisujemy wall/user/sys i jawnie oznaczamy pamięć jako `UNCONFIRMED` zamiast zgadywać.

## Poziomy akceptacji

| Poziom | Zakres | Cel |
| --- | --- | --- |
| Smoke | Mały build `--preview-rows` na jednym chromosomie | Potwierdza, że generowanie cache, lookup Fjall i podpięcie annotacji nie crashują. |
| Functional | Większy preview ograniczony do chromosomu, na przykład `chr1 preview_rows=1_000_000` | Potwierdza realistyczny payload pluginu i daje użyteczne liczby wydajnościowe. |
| Full | Pełne wejście źródłowe pluginu | Potwierdza produkcyjne generowanie cache i finalne claimy zgodności. |

## Otwarte ustalenie

Dokładna komenda `vepyr-diffly` do porównania pluginów nadal wymaga zaprojektowania. Powinna uruchamiać VEP i `vepyr` na tym samym przygotowanym VCF, włączać po jednym pluginie naraz, a potem używać istniejącego semantycznego pipeline'u porównania dla wygenerowanych annotowanych VCF-ów.
