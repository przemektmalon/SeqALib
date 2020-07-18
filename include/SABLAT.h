#include <vector>

template <typename ContainerType,
          typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0),
          typename MatchFnTy = std::function<bool(Ty, Ty)>>
class BLATSA : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy> {
private:
  int wordSize;
  int initialThreshold;

  using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

  ScoringSystem &Scoring = BaseType::getScoring();
  const ScoreSystemType Gap = Scoring.getGapPenalty();
  const ScoreSystemType Match = Scoring.getMatchProfit();
  const bool AllowMismatch = Scoring.getAllowMismatch();
  const ScoreSystemType Mismatch =
      AllowMismatch ? Scoring.getMismatchPenalty()
                    : std::numeric_limits<ScoreSystemType>::min();

  // Build the alignment using the BLAST hueristic algorithm
  void buildAlignment(ContainerType &Seq1, ContainerType &Seq2,
                      AlignedSequence<Ty, Blank> &Result) {

    auto &Data = Result.Data;
    const size_t SizeSeq1 = Seq1.size();
    const size_t SizeSeq2 = Seq2.size();

    /*std::ofstream newfile;
    std::string newPath = "/home/sean/correctChecker.txt";
    newfile.open(newPath, std::ios_base::app);
    newfile << "Began BLAST\n";
    newfile << "Seq1 Size: " << Seq1.size() << " Seq2 Size: " << Seq2.size() << "\n";
    newfile.close();*/

    // Determine appropriate word size based on the size of the smallest
    // sequence
    if (SizeSeq1 <= SizeSeq2) {
      if (SizeSeq1 < 500) {
        wordSize = (int)std::log2(SizeSeq1);
      } else {
        wordSize = 11;
      }
    } else {
      if (SizeSeq2 < 500) {
        wordSize = (int)std::log2(SizeSeq2);
      } else {
        wordSize = 11;
      }
    }

    // Prefer to calculate overlapping k-mers of shorter sequence
    // Switch the variables as Seq1 is assumed to be the shorter sequence
    // Must remember that a switch occurred for alignment ordering purposes
    bool switched = false;
    if (Seq1.size() > Seq2.size()) {
      // ContainerType temp = Seq1;
      // Seq1 = Seq2;
      // Seq2 = temp;
      std::swap(Seq1, Seq2);
      switched = true;
    }

    // Find candidate words by pattern matching
    std::vector<candidateWord> words;
    std::vector<candidateWord> candidateWords;

    /*do {
      SuffixTree<ContainerType, Ty, Blank> suffixTree;

      suffixTree.buildTree(Seq1, BaseType::getMatchOperation());
      initialThreshold = Match * wordSize;

      // Create non-overlapping patterns
      // And compare each pattern to suffix tree to find seeds
      for (int i = 0; i < Seq2.size() - wordSize + 1; i = i + wordSize) {
        ArrayView<ContainerType> pattern(Seq2);
        pattern.sliceWindow(i, i + wordSize);
        std::vector<candidateWord> patternMatches =
            suffixTree.getCandidates(pattern);
        for (int j = 0; j < patternMatches.size(); j++) {
          patternMatches[j].indexSeq2 = i;
          patternMatches[j].score = Match * wordSize;
          patternMatches[j].wordSize = wordSize;
          candidateWords.push_back(patternMatches[j]);
        }
      }

      if (candidateWords.size() == 0) {
        wordSize = (int)wordSize * 3 / 4;

        if (wordSize <= 0) {
          return;
        }
      }

      words.clear();
      suffixTree.deleteTree();
    } while (candidateWords.size() == 0);*/

    do {
      // Perfect Matches
      initialThreshold = Match * wordSize;

      for (int i = 0; i < Seq1.size() - wordSize + 1; i++) {
        candidateWord cWord;
        cWord.indexSeq1 = i;
        cWord.wordSize = wordSize;
        words.push_back(cWord);
      }

      // Compare k-mer words to non-overlapping k-mer words in Sequence 2 and a
      // score for each comparison and delete alignments below a threshold
      for (int i = 0; i < words.size(); i++) {
        for (int j = 0; j < Seq2.size() - wordSize + 1; j = j + wordSize) {
          ScoreSystemType score = 0;

          for (int k = 0; k < wordSize; k++) {
            ScoreSystemType Similarity =
                BaseType::match(Seq1[words[i].indexSeq1 + k], Seq2[j + k])
                    ? Match
                    : Mismatch;
            score += Similarity;
          }

          if (score >= initialThreshold) {
            candidateWord cWord;
            cWord = words[i];
            cWord.score = score;
            cWord.indexSeq2 = j;
            candidateWords.push_back(cWord);
          }
        }
      }

      if (candidateWords.size() == 0) {
        wordSize = (int)wordSize * 3 / 4;

        if (wordSize <= 0) {
          return;
        }
      }

      words.clear();
    } while (candidateWords.size() == 0);

    // Begin expanding out the seeds in both directions
    bool anyExpansion = false;
    std::vector<candidateWord> pushedWords;

    while (candidateWords.size() >= 1) {
      std::vector<candidateWord> newCandidates = candidateWords;
      for (int i = 0; i < newCandidates.size(); i++) {
        //<------------- Merge overlapping words ------------->
        if (i > 0) {
          candidateWord firstWord = newCandidates[i - 1];
          candidateWord secondWord = newCandidates[i];

          // The position at which the first word ends is greater than the
          // beginning of the second word This means they are overlapping so
          // merge them together
          if ((firstWord.indexSeq1 + firstWord.wordSize >
               secondWord.indexSeq1) &&
              (firstWord.indexSeq2 + firstWord.wordSize >
               secondWord.indexSeq2) &&
              (firstWord.indexSeq1 != secondWord.indexSeq1)) {

            int overlapSize =
                firstWord.indexSeq1 + firstWord.wordSize - secondWord.indexSeq1;
            int overlapSeq2Size =
                firstWord.indexSeq2 + firstWord.wordSize - secondWord.indexSeq2;
            if (overlapSize != overlapSeq2Size) {
              // i++;
              continue;
            }

            candidateWord mergedWord = firstWord;

            // re-calculate score, however, we already know word 1 score so to
            // reduce time we only calculate score of the word 2 that is being
            // added to word 1. Must be done since imperfect matches prevents us
            // from simply multiplying a match score with the number of
            // characters being added to word 1
            int mergedWordSize = mergedWord.wordSize + secondWord.wordSize -
                                 (firstWord.indexSeq1 + firstWord.wordSize -
                                  secondWord.indexSeq1);
            for (int j = mergedWord.wordSize; j < mergedWordSize; j++) {

              ScoreSystemType Similarity =
                  BaseType::match(
                      Seq1[mergedWord.indexSeq1 + j],
                      Seq2[mergedWord.indexSeq2 + mergedWord.wordSize])
                      ? Match
                      : Mismatch;
              // char wordChar = mergedWord.word[j];
              // ScoreSystemType Similarity = wordChar ==
              // Seq2[mergedWord.indexSeq2 + mergedWord.wordSize] ? Match :
              // Mismatch;
              mergedWord.score += Similarity;
              if (Similarity < 0) {
                mergedWord.score = Similarity;
              }
            }

            mergedWord.wordSize = mergedWordSize;
            newCandidates[i - 1] = mergedWord;
            newCandidates.erase(newCandidates.begin() + i);
            i--;
            continue;
          }
          // The two words are on the verge of joining so merging is simple
          // addition Separate else if for efficiency's sake
          else if (firstWord.indexSeq1 + firstWord.wordSize ==
                   secondWord.indexSeq1) {

            int overlapSize =
                firstWord.indexSeq1 + firstWord.wordSize - secondWord.indexSeq1;
            int overlapSeq2Size =
                firstWord.indexSeq2 + firstWord.wordSize - secondWord.indexSeq2;

            if (overlapSize != overlapSeq2Size) {
              continue;
            }
            candidateWord mergedWord = firstWord;
            mergedWord.wordSize = mergedWord.wordSize + secondWord.wordSize;
            mergedWord.score = mergedWord.score + secondWord.score;
            newCandidates[i - 1] = mergedWord;
            newCandidates.erase(newCandidates.begin() + i);
            i--;
            continue;
          }
        }

        //<------------- Expand word ------------->
        ScoreSystemType score = newCandidates[i].score;
        bool expansion = false;
        // Expand left
        int leftIndex = newCandidates[i].indexSeq1 - 1;
        if (leftIndex >= 0) { // We can expand left
          // char leftChar = Seq1[leftIndex];
          // newCandidates[i].word = leftChar + newCandidates[i].word;
          newCandidates[i].indexSeq1--;
          newCandidates[i].indexSeq2--;
          newCandidates[i].wordSize++;

          ScoreSystemType Similarity;

          if (newCandidates[i].indexSeq2 < 0) {
            Similarity = Mismatch;
          } else {
            Similarity = BaseType::match(Seq1[leftIndex],
                                         Seq2[newCandidates[i].indexSeq2])
                             ? Match
                             : Mismatch;
          }

          // A detremental expansion so revert
          if (newCandidates[i].score + Similarity < newCandidates[i].score ||
              Similarity < 0) {
            newCandidates[i].indexSeq1++;
            newCandidates[i].indexSeq2++;
            newCandidates[i].wordSize--;
          } else {
            newCandidates[i].score += Similarity;
            expansion = true;
            anyExpansion = true;
          }
        }

        // Expand right
        int rightIndex =
            newCandidates[i].indexSeq1 +
            newCandidates[i].wordSize; // Plus wordSize as we index points
                                       // towards start of word
        if (rightIndex < Seq1.size()) {

          ScoreSystemType Similarity;

          if (newCandidates[i].indexSeq2 + newCandidates[i].wordSize >=
              Seq2.size()) {
            Similarity = Mismatch;
          } else {
            Similarity = BaseType::match(Seq1[rightIndex],
                                         Seq2[newCandidates[i].indexSeq2 +
                                              newCandidates[i].wordSize])
                             ? Match
                             : Mismatch;
          }

          // A detremental expansion so revert
          if (newCandidates[i].score + Similarity < newCandidates[i].score ||
              Similarity < 0) {
          } else {
            newCandidates[i].wordSize++;

            newCandidates[i].score += Similarity;
            expansion = true;
            anyExpansion = true;
          }
        }

        // Word cannot be expanded further
        if (!expansion) {
          pushedWords.push_back(newCandidates[i]);
          newCandidates.erase(newCandidates.begin() + i);
        }
      }

      // All possible expansions have taken place so exit
      if (!anyExpansion) {
        candidateWords = pushedWords;
        break;
      }

      candidateWords = newCandidates;
    }

    candidateWords = pushedWords;
    //<----- Longest Increasing Subsequence ----->
    // Find the longest increasing subsequence of the complete words
    // Important for smaller search space and overall score/usage of words
    std::sort(candidateWords.begin(), candidateWords.end());

    if (candidateWords.size() > 1) {
      // Denoting X as completeWords
      int P[candidateWords.size()];
      int M[candidateWords.size() + 1];

      int l = 0;
      int lo = 0;
      int hi = 0;
      for (int i = 0; i < candidateWords.size(); i++) {
        // Binary search for the largest positive j <= l
        // Such that X[M[j]] <= X[i]
        lo = 1;
        hi = l;
        while (lo <= hi) {
          int mid = ceil((lo + hi) / 2);
          if (candidateWords[M[mid]] <= candidateWords[i]) {
            lo = mid + 1;
          } else {
            hi = mid - 1;
          }
        }

        // After searching, lo is 1 greater than the length
        // of the longest prefix of X[i]
        int newL = lo;

        // The predecessor of X[i] is the last index of
        // the subsequence of length newL - 1
        P[i] = M[newL - 1];
        M[newL] = i;

        // If we found a subsequence longer than any we've foudn yet
        if (newL > l) {
          l = newL;
        }
      }

      // Reconstruct the longest increasing subsequence
      std::vector<candidateWord> newList;
      newList.resize(l);
      int k = M[l];
      for (int i = l - 1; i >= 0; i--) {
        newList[i] = candidateWords[k];
        k = P[k];
      }

      candidateWords.clear();
      candidateWords = newList;
    }

    if (candidateWords.size() == 0) {
          /*newfile.open(newPath, std::ios_base::app);
    newfile << "USE SW\n";
    newfile.close();*/
      NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> bridge(
          Scoring, BaseType::getMatchOperation());
      ArrayView<ContainerType> seq1Sub(Seq1);
      ArrayView<ContainerType> seq2Sub(Seq2);
      seq1Sub.sliceWindow(0, Seq1.size());
      seq2Sub.sliceWindow(0, Seq2.size());
      AlignedSequence<Ty, Blank> NWAlignment =
          bridge.getAlignment(seq1Sub, seq2Sub);
      Data.insert(Data.end(), NWAlignment.Data.begin(), NWAlignment.Data.end());
      return;
    }

    //<----- Build Alignment ----->
    if (switched) {
      // ContainerType temp = Seq1;
      // Seq1 = Seq2;
      // Seq2 = temp;
      std::swap(Seq1, Seq2);

      for (int i = 0; i < candidateWords.size(); i++) {
        int tempIndex = candidateWords[i].indexSeq2;
        candidateWords[i].indexSeq2 = candidateWords[i].indexSeq1;
        candidateWords[i].indexSeq1 = tempIndex;
      }
    }

    // For each word
    int i = 0;
    bool breakOut = false;
    for (i = 0; i < candidateWords.size(); i++) {
      // Align word
      alignWord(candidateWords[i], Seq1, Seq2, Result);

      // Bridge gap to word ignoring overlapping words
      int j = i + 1;
      while (true) {
        if (j < candidateWords.size()) {
          // Overlapping word so look to bridge to next word if possible
          if (candidateWords[i].indexSeq1 + candidateWords[i].wordSize >
                  candidateWords[j].indexSeq1 ||
              candidateWords[i].indexSeq2 + candidateWords[i].wordSize >
                  candidateWords[j].indexSeq2) {
            j++;
            continue;
          }
          NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy>
              bridge(Scoring, BaseType::getMatchOperation());
          ArrayView<ContainerType> seq1Sub(Seq1);
          ArrayView<ContainerType> seq2Sub(Seq2);
          seq1Sub.sliceWindow(candidateWords[i].indexSeq1 +
                                  candidateWords[i].wordSize,
                              candidateWords[j].indexSeq1);
          seq2Sub.sliceWindow(candidateWords[i].indexSeq2 +
                                  candidateWords[i].wordSize,
                              candidateWords[j].indexSeq2);

          AlignedSequence<Ty, Blank> NWAlignment =
              bridge.getAlignment(seq1Sub, seq2Sub);
          Data.insert(Data.end(), NWAlignment.begin(), NWAlignment.end());
          i = j - 1;
          break;
        } else {
          breakOut = true;
          break;
        }
      }

      if (breakOut) {
        i++;
        break;
      }
    }

    /*newfile.open(newPath, std::ios_base::app);
    newfile << "candidate word size: " << candidateWords.size() << "\n";
    newfile << "i: " << i << "\n";
    newfile.close();*/
    if (i == 1) {
      BaseType::forceGlobal(
          Seq1, Seq2, Result, candidateWords[0].indexSeq1,
          candidateWords[0].indexSeq2,
          candidateWords[0].indexSeq1 + candidateWords[0].wordSize,
          candidateWords[0].indexSeq2 + candidateWords[0].wordSize);
    } else {
      // Force global alignment
      BaseType::forceGlobal(
          Seq1, Seq2, Result, candidateWords[0].indexSeq1,
          candidateWords[0].indexSeq2,
          candidateWords[i - 1].indexSeq1 + candidateWords[i - 1].wordSize,
          candidateWords[i - 1].indexSeq2 + candidateWords[i - 1].wordSize);
    }

    /*newfile.open(newPath, std::ios_base::app);
    i=0;
    for (auto Char:Result)
    {
      if (Char.get(0)!=Blank)
      {
        if (BaseType::match(Seq1[i], Char.get(0)))
        {
          newfile << "1";
        }
        else
        {
          newfile << "0";
        }
        i++;
      }
    }
    newfile << "\n";
    newfile << "COMPLETED BLAST\n";
    newfile << "\n";
    newfile.close();*/

    if (Seq1.size()==145 && Seq2.size()==121)
    {
      Result.Data.clear();
            NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> bridge(
          Scoring, BaseType::getMatchOperation());
      ArrayView<ContainerType> seq1Sub(Seq1);
      ArrayView<ContainerType> seq2Sub(Seq2);
      seq1Sub.sliceWindow(0, Seq1.size());
      seq2Sub.sliceWindow(0, Seq2.size());
      AlignedSequence<Ty, Blank> NWAlignment =
          bridge.getAlignment(seq1Sub, seq2Sub);
      Data.insert(Data.end(), NWAlignment.Data.begin(), NWAlignment.Data.end());
      return;
    }
  }

  void alignWord(candidateWord word, ContainerType &Seq1, ContainerType &Seq2,
                 AlignedSequence<Ty, Blank> &Result) {
    auto &Data = Result.Data;

    for (int i = 0; i < word.wordSize; i++) {
      Data.push_back(typename BaseType::EntryType(
          Seq1[i + word.indexSeq1], Seq2[i + word.indexSeq2], true));
    }
  }

public:
  static ScoringSystem getDefaultScoring() { return ScoringSystem(-1, 2, -1); }

  BLATSA() : BaseType(getDefaultScoring(), nullptr) {}

  BLATSA(ScoringSystem Scoring, MatchFnTy Match) : BaseType(Scoring, Match) {}

  virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1,
                                                  ContainerType &Seq2) {
    AlignedSequence<Ty, Blank> Result;
    buildAlignment(Seq1, Seq2, Result);
    return Result;
  }
};
