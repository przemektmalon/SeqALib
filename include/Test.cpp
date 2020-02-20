#include <iostream>
#include <string>
#include <chrono>

#include "SequenceAlignment.h"

template <typename T>
bool equal(T V1, T V2) { return V1 == V2; }

template <typename Ty, Ty Blank>
void printAlignment(AlignedSequence<Ty, Blank> &Result)
{
  for (auto &Entry : Result)
  {
    std::cout << Entry.get(0);
  }
  std::cout << std::endl;
  for (auto &Entry : Result)
  {
    if (Entry.match())
      std::cout << "|";
    else
      std::cout << " ";
  }
  std::cout << std::endl;
  for (auto &Entry : Result)
  {
    std::cout << Entry.get(1);
  }
  std::cout << std::endl;
}

int main(int argc, char **argv)
{

  //std::string seq1 = "AATCG";
  //std::string seq2 = "AACG";
  //std::string seq1 = "AGGATCGGCTAGAGCTAGAGCTAGCTAGTAGC";
  //std::string seq2 = "GAGATCGGCGGATTACAGGCTATCGA";
  //std::string seq1 = "AAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
  //std::string seq2 = "AAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
  //std::string seq1 = "AAAAGGGGTTTTCCCC";
  //std::string seq2 = "AAAAGGGGTTTTCCCC";
  //std::string seq1 = "AAAA";
  //std::string seq2 = "AAAA";
  std::string seq1 = "CTGAAGCGG";
  std::string seq2 = "CTCAAGCGTAGTCC";
  //std::string seq1 = "AGTAC";
  //std::string seq2 = "AAG";
  //std::string seq1 = "GGTCGCGACTACGTGAGCTAGGGCTCCGGACTGGGCTGTATAGTCGAGTC";
  //std::string seq2 = "TGATCTCGCCCCGACAACTGCAAACCCCAACTTATTTAGATAACATGGTTTACAGG";
  //std::string seq1 = "AGGGGGGGGTTTTCAAAAGCTCTTCGCATGCGCATAGGCTTAGGAGTCGAGGATCGGCGGCATATTAAGAGAGGGCGGCGCCCATATATTAGGCGCGCGCGATTATATTATAATATATATATTAATCGGGCGCTATGCGC";
  //std::string seq2 = "AGGTCTCCGGCATGAGCTAGCGGCGAGTTATAGCGCGTTTCAAAAGCTCTTCGCATGCGCATAGGCTTAGGAGTCGAGGATCGGCGGCATATTAAGAGAGGGCGGCAGGCTTGGAAAAAGGCTCTGAGATGCGAGAATAGAGGAGAG";
  //std::string seq1 = "GTAGATTGCCATGCGTAGAGCTAACGAGCCAGCGGAAAGCGTGAGGCGCTTTTAAGCATGGCGAGTAAGTGATCCAACGCTTCGGATATGACTATATACTTAGGTTCGATCTCGTCCCGAGAATTCTAAGCCTCAACATCTATGAGTTATGAGGTTAGCCGAAAAAGCACGTGGTGGCGCCCACCGACTGTTCCCAGACTGTAGCTCTTTGTTCTGTCAAGGCCCGACCTTCATCGCGGCCGATTCCTTCTGCGGACCATACCGTCCTGATACTTTGGTCATGTTTCCGTTGTAGGAGTGAACCCACTTGCCTTTGCGTCTTAATACCAATGAAAAACCTATGCACTTTGTACAGGGTACCATCGGGATTCTGAACCCTCAGATAGTGGGGATCCCGGGTATAGACCTTTATCTGCGGTCCAACTTAGGCATAAACCTCCATGCTACCTTGTCAGACCCACCCTGCACGAGGTAAATATGGGACGCGTCCGACCTGGCTCCTGGCGTTCTACGCCGCCACGTGTTCGTTAACTGTTGATTGGTAGCACAAAAGTAATACCATGGTCCTTGAAATTCGGCTCAGTTAGTTCGAGCGTAATGTCACAAATGGCGCAGAACGGCAATGAGTGTTTGACACTAGGTGGTGTTCAGTTCGGTAACGGAGAGACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCTAGGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTCGATTGAAGATGGGTAGGGGGCTTCAATTCGTCCAGCACTCTACGGTACCTCCGAGAGCAAGTAGGGCACCCTGTAGTTCGAAGCGGAACTATTTCGTGGGGCGAGCCCACATCGTCTCTTCTGCGGATGACTTAACACGTTAGGGAGGTGGAGTTGATTCGAACGATGGTTATAAATCAAAAAAACGGAACGCTGTCTGGAGGATGAATCTAACGGTGCGTAACTCGATCACTCGTAGATTGCCATGCGTAGAGCTAACGAGCCAGCGGAAAGCGTGAGGCGCTTTTAAGCATGGCGAGTAAGTGATCCAACGCTTCGGATATGACTATATACTTAGGTTCGATCTCGTCCCGAGAATTCTAAGCCTCAACATCTATGAGTTATGAGGTTAGCCGAAAAAGCACGTGGTGGCGCCCACCGACTGTTCCCAGACTGTAGCTCTTTGTTCTGTCAAGGCCCGACCTTCATCGCGGCCGATTCCTTCTGCGGACCATACCGTCCTGATACTTTGGTCATGTTTCCGTTGTAGGAGTGAACCCACTTGCCTTTGCGTCTTAATACCAATGAAAAACCTATGCACTTTGTACAGGGTACCATCGGGATTCTGAACCCTCAGATAGTGGGGATCCCGGGTATAGACCTTTATCTGCGGTCCAACTTAGGCATAAACCTCCATGCTACCTTGTCAGACCCACCCTGCACGAGGTAAATATGGGACGCGTCCGACCTGGCTCCTGGCGTTCTACGCCGCCACGTGTTCGTTAACTGTTGATTGGTAGCACAAAAGTAATACCATGGTCCTTGAAATTCGGCTCAGTTAGTTCGAGCGTAATGTCACAAATGGCGCAGAACGGCAATGAGTGTTTGACACTAGGTGGTGTTCAGTTCGGTAACGGAGAGACTGTGCGGCATACTTAATTATACATTTGAAACGCGCCCAAGTGACGCTAGGCAAGTCAGAGCAGGTTCCCGTGTTAGCTTAAGGGTAAACATACAAGTCGATTGAAGATGGGTAGGGGGCTTCAATTCGTCCAGCACTCTACGGTACCTCCGAGAGCAAGTAGGGCACCCTGTAGTTCGAAGCGGAACTATTTCGTGGGGCGAGCCCACATCGTCTCTTCTGCGGATGACTTAACACGTTAGGGAGGTGGAGTTGATTCGAACGATGGTTATAAATCAAAAAAACGGAACGCTGTCTGGAGGATGAATCTAACGGTGCGTAACTCGATCACTC";
  //std::string seq2 = "ACTCGCTATTCGAACTGCGCGAAAGTTCCCAGCGCTCATACACTTGGTTCCGAGGCCTGTCCTGATATATGAACCCAAACTAGAGCGGGGCTGTTGACGTTTGGAGTTGAAAAAATCTAATATTCCAATCGGCTTCAACGTGCACCACCGCAGGCGGCTGACGAGGGGCTCACACCGAGAAAGTAGACTGTTGCGCGTTGGGGGTAGCGCCGGCTAACAAAGACGCCTGGTACAGCAGGAGTATCAAACCCGTACAAAGGCAACATCCTCACTTCGGTGAATCGAAGCGCGGCATCAGGGTTACTTTTTGGATACCTGAAACAAAACCCATCGTAGTCCTTAGACTTGGCACACTTACACCTGCAGCGCGCGCATCTGGAAATAGAGGCCAAGTTCGATCCGTACTCCGACGTACGATGCAACAGTGTGGATGTGACGAGCTTCATTTATACCCTTCGCGCGCCGGACCGGCCTCCGCAAGGCGCGGCGGTGCACAAGCAATTGACAACTAACCACCGTGTATTCGTTATGGCATCAGGCAGTTTAAGTCGAGACAATAGGGCTCGCAATACACAGTTTACCGCATCTTGCCCTAACTGACAAACTGTGATCGACCACTAGCCATGCCATTGCCTCTTAGACACCCCGATACAGTGATTATGAAAGGTTTGCGGGGCATGGCTACGACTTGTTCAGCTACGTCCGAGGGCAGAAACTTATCCCCATTTGTATGTTCACCTATCTACTACCCATCCCCGGAGGTTAAGTAGGTTGTGAGATGCGGGAGAGGTTCTCGATCTTCCCGTGGGACGTCAACCTTTCCCTTGATAAAGCATCCCGCTCGGGTATGGCAGTGAGTACGCCTTCTGAATTGTGCTATCCTTCGTCCTTATCAAAGCTTGCTACCAATAATTAGGATTATTGCCTTGCGACAGACTTCCTACTCACACTCCCTCACATTGAGCTACTCGATGGGCGATTAGCTTGACCCGCTCTGTAGACTCGCTATTCGAACTGCGCGAAAGTTCCCAGCGCTCATACACTTGGTTCCGAGGCCTGTCCTGATATATGAACCCAAACTAGAGCGGGGCTGTTGACGTTTGGAGTTGAAAAAATCTAATATTCCAATCGGCTTCAACGTGCACCACCGCAGGCGGCTGACGAGGGGCTCACACCGAGAAAGTAGACTGTTGCGCGTTGGGGGTAGCGCCGGCTAACAAAGACGCCTGGTACAGCAGGAGTATCAAACCCGTACAAAGGCAACATCCTCACTTCGGTGAATCGAAGCGCGGCATCAGGGTTACTTTTTGGATACCTGAAACAAAACCCATCGTAGTCCTTAGACTTGGCACACTTACACCTGCAGCGCGCGCATCTGGAAATAGAGGCCAAGTTCGATCCGTACTCCGACGTACGATGCAACAGTGTGGATGTGACGAGCTTCATTTATACCCTTCGCGCGCCGGACCGGCCTCCGCAAGGCGCGGCGGTGCACAAGCAATTGACAACTAACCACCGTGTATTCGTTATGGCATCAGGCAGTTTAAGTCGAGACAATAGGGCTCGCAATACACAGTTTACCGCATCTTGCCCTAACTGACAAACTGTGATCGACCACTAGCCATGCCATTGCCTCTTAGACACCCCGATACAGTGATTATGAAAGGTTTGCGGGGCATGGCTACGACTTGTTCAGCTACGTCCGAGGGCAGAAACTTATCCCCATTTGTATGTTCACCTATCTACTACCCATCCCCGGAGGTTAAGTAGGTTGTGAGATGCGGGAGAGGTTCTCGATCTTCCCGTGGGACGTCAACCTTTCCCTTGATAAAGCATCCCGCTCGGGTATGGCAGTGAGTACGCCTTCTGAATTGTGCTATCCTTCGTCCTTATCAAAGCTTGCTACCAATAATTAGGATTATTGCCTTGCGACAGACTTCCTACTCACACTCCCTCACATTGAGCTACTCGATGGGCGATTAGCTTGACCCGCTCTGTAG";
  //std::string seq1 = "ACGGTTGC";
  //std::string seq2 = "AGCGTC";
  //std::string seq1 = "AGTC";
  //std::string seq2 = "GACT";
  //std::string seq1 = "AAGG";
  //std::string seq2 = "AAGG";
  //std::string seq1 = "AAGGGGGCAACCCAATTGTCAAAA";
  //std::string seq2 = "AAGTTGGCGGCCCAAGCTGCGAAA";
  //std::string seq1 = "AAAGGGTTTCCCAAGGTTCCAGTC";
  //std::string seq2 = "AAAGGGTTTCCCAAGGCTCCAGTC";

  // AAA GAA TGCAT
  // | |  |  | |||
  // A A  A CT CAT

  // AAA GAATGCAT
  // |||    | |||
  // AAAC   T CAT

  if (argc > 1)
  {
    seq1 = std::string(argv[1]);
    seq2 = std::string(argv[2]);
  }

  ScoringSystem scoringSystem = ScoringSystem(-2, 1, -1, false);
  ScoringSystem affineSystem = ScoringSystem(-3, -1, 1, -1, false);

  auto t1 = std::chrono::high_resolution_clock::now();
  NeedlemanWunschSA<std::string, char, '-'> SA(scoringSystem, equal<char>); //could also use `nullptr` instead of `equal<char>`
  AlignedSequence<char, '-'> NWAlignment = SA.getAlignment(seq1, seq2);
  auto t2 = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

  std::cout << "# Needleman-Wunsch:" << std::endl;
  printAlignment(NWAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;

  /*auto t3 = std::chrono::high_resolution_clock::now();
  AlignedSequence<char,'-'> HAlignment = HirschbergSA<std::string,char,'-'>(scoringSystem,equal<char>).getAlignment(seq1,seq2);
  auto t4 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
  std::cout << "# Hirschberg:" << std::endl;
  printAlignment(HAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;*/

  auto t5 = std::chrono::high_resolution_clock::now();
  AlignedSequence<char, '-'> SWAlignment = SmithWatermanSA<std::string, char, '-'>(scoringSystem, equal<char>).getAlignment(seq1, seq2);
  auto t6 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();
  std::cout << "# Smith-Waterman:" << std::endl;
  printAlignment(SWAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;

  auto t7 = std::chrono::high_resolution_clock::now();
  AlignedSequence<char, '-'> LGAlignment = LocalGotohSA<std::string, char, '-'>(affineSystem, equal<char>).getAlignment(seq1, seq2);
  auto t8 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();
  std::cout << "# Local Gotoh:" << std::endl;
  printAlignment(LGAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;

  auto t9 = std::chrono::high_resolution_clock::now();
  AlignedSequence<char, '-'> GGAlignment = GlobalGotohSA<std::string, char, '-'>(affineSystem, equal<char>).getAlignment(seq1, seq2);
  auto t10 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(t10 - t9).count();
  std::cout << "# Global Gotoh:" << std::endl;
  printAlignment(GGAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;

  auto t11 = std::chrono::high_resolution_clock::now();
  AlignedSequence<char, '-'> MMAlignment = MyersMillerSA<std::string, char, '-'>(affineSystem, equal<char>).getAlignment(seq1, seq2);
  auto t12 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(t12 - t11).count();
  std::cout << "# Myers Miller: " << std::endl;
  printAlignment(MMAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;

  /*auto t15 = std::chrono::high_resolution_clock::now();
  AlignedSequence<char, '-'> FOGSAAAlignment = FOGSAASA<std::string, char, '-'>(scoringSystem, equal<char>).getAlignment(seq1, seq2);
  auto t16 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(t16 - t15).count();
  std::cout << "# FOGSAA: " << std::endl;
  printAlignment(FOGSAAAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;

  auto t17 = std::chrono::high_resolution_clock::now();
  AlignedSequence<char, '-'> BLATAlignment = BLATSA<std::string, char, '-'>(scoringSystem, equal<char>).getAlignment(seq1, seq2);
  auto t18 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(t18 - t17).count();
  std::cout << "# BLAT: " << std::endl;
  printAlignment(BLATAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;

  auto t19 = std::chrono::high_resolution_clock::now();
  AlignedSequence<char, '-'> GappedBLATAlignment = GappedBLATSA<std::string, char, '-'>(scoringSystem, equal<char>).getAlignment(seq1, seq2);
  auto t20 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(t20 - t19).count();
  std::cout << "# Gapped BLAT: " << std::endl;
  printAlignment(GappedBLATAlignment);
  std::cout << "Time Taken: " << duration << " microseconds" << std::endl;
  std::cout << "" << std::endl;*/

  /*SuffixTree<std::string, char, '-'> suffixTree;
  std::string seq = seq1 + "#" + seq2;
  suffixTree.buildTree(seq, seq1.size() + 1, SA.getMatchOperation());
  std::cout << "boutta print" << std::endl;
  suffixTree.getLCS();*/

  //std::getchar();

  return 0;
}
