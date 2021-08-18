package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

var hg38size int = 3272116950

const (
	krit1 = "seqs/krit1_exons.fasta"
	ccm2 = "seqs/ccm2_exons.fasta"
	pdcd10 = "seqs/pdcd10_exons.fasta"
)

func numLofSnv(filename string) int {
	f := fasta.Read(filename)
	exonIntronJunctions := (len(f) * 2) - 2
	var seq []dna.Base
	for _, fseq := range f {
		seq = append(seq, fseq.Seq...)
	}

	var possibleLofSnv int = exonIntronJunctions * 4 * 3 // 4 for +2/-2 canonical sites, 3 for all base but ref

	codons := dna.BasesToCodons(seq)
	for _, c := range codons {
		possibleLofSnv += lofSnvPerCodon(c)
	}

	return possibleLofSnv
}

func lofSnvPerCodon(c dna.Codon) int {
	var possibleLofSnv int
	if dna.GeneticCode[c] == dna.Stop {
		return 0
	}

	testCodons := allSnvCodons(c)
	for _, tc := range testCodons {
		if dna.GeneticCode[tc] == dna.Stop {
			possibleLofSnv++
		}
	}

	return possibleLofSnv
}

func allSnvCodons(c dna.Codon) []dna.Codon {
	var answer []dna.Codon
	for i := 0; i < 3; i++ {
		orig := c
		c[i] = dna.A
		if c != orig {
			answer = append(answer, c)
		}
		c[i] = dna.C
		if c != orig {
			answer = append(answer, c)
		}
		c[i] = dna.G
		if c != orig {
			answer = append(answer, c)
		}
		c[i] = dna.T
		if c != orig {
			answer = append(answer, c)
		}
		c = orig
	}
	return answer
}

func main() {
	fmt.Println("KRIT1  LOF: ", numLofSnv(krit1))
	fmt.Println("CCM2   LOF: ", numLofSnv(ccm2))
	fmt.Println("PDCD10 LOF: ", numLofSnv(pdcd10))
}
