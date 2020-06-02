package main

import (
  "fmt"
  "io"
  "log"
  "os"
  "strings"

  "github.com/biogo/hts/bam"
  "github.com/biogo/hts/sam"
)

func DecodeQual(qual []byte) string {
  squal := make([]string, 0, len(qual))
  for pos := range qual {
    squal = append(squal, string(qual[pos]+33))
  }
  return strings.Join(squal, "")
}

func main() {
  bam_file := os.Args[1]

  f, err := os.Open(bam_file)
  if err != nil {
    log.Fatalf("could not open file %q:", err)
  }
  defer f.Close()

  b, err := bam.NewReader(f, 1)
  if err != nil {
    log.Fatalf("could not read bam:", err)
  }
  defer b.Close()

  var read_name string
  for {
    rec, err := b.Read()
    if err == io.EOF {
      break
    }
    if err != nil {
      log.Fatalf("error reading bam: %v", err)
    }
    if rec != nil {
      var read []string

      // See if the cell barcode tag is available
      cb := rec.AuxFields.Get(sam.NewTag("CB"))
      if cb != nil {
        read_name = strings.Join([]string{"@", rec.Name, "_", cb.String()}, "")
      } else {
        read_name = strings.Join([]string{"@", rec.Name}, "")
      }

      // See if the molecular barcode is available
      ub := rec.AuxFields.Get(sam.NewTag("UB"))
      if ub != nil {
        read_name = strings.Join([]string{read_name, "_", ub.String()}, "")
      }

      read = append(read, read_name)
      read = append(read, string(rec.Seq.Expand()))
      read = append(read, "+")
      read = append(read, DecodeQual(rec.Qual))

      fmt.Println(strings.Join(read, "\n"))
    }
  }

}
