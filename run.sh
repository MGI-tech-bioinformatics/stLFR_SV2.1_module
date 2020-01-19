#!/bin/bash
perl ./LFR-sv.pl -bam example/L0.sort.rmdup.bam -out example/result -ncpu 20 -phase example/phase_out -bl data/bl_region_hg19_nochr -cl data/con_hg19_nochr -human Y
