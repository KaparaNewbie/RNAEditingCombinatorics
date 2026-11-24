# Historical software for accepted manuscript stages

Use the [repository-root README](../README.md) for maintained setup and the [software directory guide](README.md) for current specifications, dated inventories and [archived historical environments](old_versions/README.md). The earlier development commands are preserved in [Logs.Final.md](../Logs.Final.md). Archiving those specifications does not establish their exact use in any stage below.

Current installations and historical execution are separate evidence. Native headers were read directly on 23 November 2026, without traversing BAM reads or appending inspection program records. The review archive contains `native_headers.json`, `july_native_provenance.json`, version-command receipts and the complete current inventories.

| Stage | Historical version(s) | Evidence and limits |
|---|---|---|
| Original Illumina alignment/filtering | bwa-mem2 2.2.1; samtools 1.14 | Native `D.pealeii/Alignment/Illumina/reads.sorted.aligned.filtered.bam` program records. A later `samtools view -H` command in old audit text is inspection, not the processing version. |
| PB1 consensus, alignment and quality filtering | CCS 6.2.0; pbmm2 1.7.0; samtools 1.18 | Native `D.pealeii/MpileupAndTranscripts/RQ998.TopNoisyPositions3.BQ30/PCLO-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.bam` header; original stage uses the original squid reference. |
| July 2026 PB1 and Illumina complete-noise processing | samtools 1.18 / HTSlib 1.17 | Native July 21 execution logs in `D.pealeii/MpileupAndTranscripts/CompleteNoisePositions/` name the `combinatorics/bin/samtools` executable; native PB1 BAM headers record 1.18; current retained executable reports 1.18 / 1.17. Historical distinction was verified in the preceding integration. This task reread logs/headers, not the historical executable bytes. Uses Jan2025 reference/annotations. |
| PB3 delivered HiFi/CCS reads (JR1–JR3) | Provider CCS 8.2.0; lima 2.11.0; pbtrim 1.1.0; jasmine 2.3.0 | All three native delivered BAM headers in `D.pealeii/Data/AdditionalRawWithUMIs/CCS/`. These are delivered reads, not a request to rerun provider processing. |
| PB3 March 2026 alignment/merge/quality filtering | pbmm2 26.1.99; samtools 1.22.1 | Native `D.pealeii/MpileupAndTranscripts/AdditionalUMILongReads/comp134400_c0_seq1_extended.merged.MinRQ998.bam` header. Jan2025 extended reference. Exact uncommitted March source snapshot remains unrecovered. |
| Previously published O. vulgaris polished-read handoff/refinement | CCS 4.0.0; lima 2.6.0; Iso-Seq refine 3.8.2; pbmerge 1.0.6 | Native `O.vulgaris/Data/PRJNA791920/IsoSeqPolished/Refined/SRR17321895.flnc.bam` header; one sampled run, not verification of every historical run. |
| November accepted final analyses | Python 3.13.14; NumPy 2.4.6; pandas 3.0.3; scikit-learn 1.9.0 | Existing `combinatorics2` metadata plus prior integration runtime receipts; current S16 executes in this prefix. |
| Julia analyses / S19 producer | Julia 1.7.2; tracked Project.toml / Manifest.toml dependency pins | Installed executable version, read-only manifest/source inspection and previous S19 runtime receipts. A single Julia version does not imply every historical package tree is byte-identical to today's depot. |

Trimmomatic 0.39 and PRINSEQ-lite 0.20.4 remain historical environment candidates, not independently established execution versions. Instrument/basecaller software for the original Illumina delivery and full acquisition provenance require provider records. Do not replace those gaps with current installed versions.
