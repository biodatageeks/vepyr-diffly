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
| `spliceai` | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `cadd` | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | CADD ma dwa wejścia źródłowe: SNV i indel. Domyślnie zapisujemy sumaryczne rozmiary wejścia/cache, chyba że testujemy subźródła osobno. |
| `alphamissense` | Full cache parquet-only + plugin E2E on `/tmp/hg002_grch38_100k_stratified.vcf` | 661.0s convert / 661.0s wall | 591.4MB gz / 1.1GB raw | 2.0G | n/a (`--no-plugin-fjall`) | tak (plugin-only compare, 100k) | ~896s wall (`2026-04-20T08:26:17+02:00` -> `2026-04-20T08:41:13+02:00`) | ~759s wall (`2026-04-20T08:41:13+02:00` -> `2026-04-20T08:53:52+02:00`) | 2026-04-20: `runs/alphamissense-e2e-100k-bench/summary.json` dał `variant.equal=true`, `consequence.equal=true`, `427` joined-equal consequence rows, `0/0/0` mismatchy przy `--compare-only-plugins`. Outputy annotacji: `vep.annotated.vcf=579M` (`101069` linii), `vepyr.annotated.vcf=585M` (`101063` linii). Dodatni smoke fixture `/tmp/alphamissense_positive_10.clean.vcf` też dał `10/10` joined-equal. Peak RSS: `UNCONFIRMED`. |
| `dbnsfp` | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | dbNSFP prawdopodobnie będzie najdroższą ścieżką przygotowania źródła dla preview i full-cache. |

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
