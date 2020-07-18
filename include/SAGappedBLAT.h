#include <vector>
#include <queue>
template <typename ContainerType,
          typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0),
          typename MatchFnTy = std::function<bool(Ty, Ty)>>
class GappedBLATSA
    : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy> {
private:
  int wordSize;
  int initialThreshold;
  int distance;
  int gapExpansionThreshold = 0;

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
    size_t SizeSeq1 = Seq1.size();
    size_t SizeSeq2 = Seq2.size();

    // Create diagonal vector that stores index of most recent hit on that
    // diagonal
    std::vector<int> diagonalVector;
    std::vector<int>
        diagonalWordVector; // Stores word size of most recent diagonal so we
                            // can ignore overlapping words

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

      std::swap(SizeSeq1, SizeSeq2);
    }

    diagonalVector.resize(SizeSeq1 + SizeSeq2 + 1);
    diagonalWordVector.resize(SizeSeq1 + SizeSeq2 + 1);
    std::fill(diagonalVector.begin(), diagonalVector.end(), -1);
    std::fill(diagonalWordVector.begin(), diagonalWordVector.end(), -1);

    distance = 4 * wordSize;

    // Prefer to calculate overlapping k-mers of shorter sequence
    // Switch the variables as Seq1 is assumed to be the shorter sequence
    // Must remember that a switch occurred for alignment ordering purposes
    bool switched = false;
    if (SizeSeq1 > SizeSeq2) {
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
        // std::string string = Seq2.substr(i, wordSize);
        // std::cout << "PATTERN MATCHING: " << string << std::endl;
        ArrayView<ContainerType> pattern(Seq2);
        pattern.sliceWindow(i, i + wordSize);
        std::vector<candidateWord> patternMatches =
            suffixTree.getCandidates(pattern);
        for (int j = 0; j < patternMatches.size(); j++) {
          // std::string stringa = Seq1.substr(patternMatches[j].indexSeq1,
          // wordSize);

          // std::cout << "A candidate: " << stringa << std::endl;
          patternMatches[j].indexSeq2 = i;
          patternMatches[j].score = Match * wordSize;
          patternMatches[j].wordSize = wordSize;
          candidateWords.push_back(patternMatches[j]);
        }
        // std::cout << std::endl;
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
      initialThreshold = Match * wordSize;

      for (int i = 0; i < Seq1.size() - wordSize + 1; i++) {
        candidateWord cWord;
        cWord.indexSeq1 = i;
        cWord.wordSize = wordSize;
        words.push_back(cWord);
      }

      // Compare k-mer words to non-overlapping k-mer words in Sequence 2
      // and a score for each comparison and delete alignments below a threshold
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

    gapExpansionThreshold = (wordSize * Match) + (Match * log2(SizeSeq2)) - 5;

    std::vector<candidateWord> completeWords;

    // Begin expanding out the seeds in both directions
    while (candidateWords.size() >= 1) {
      std::vector<candidateWord> pushedWords;

      bool anyExpansion = false;
      for (int i = 0; i < candidateWords.size(); i++) {

        // Only expand if within distance A of the previous hit
        bool allowExtension = false;
        int diagonal = (SizeSeq2 + candidateWords[i].indexSeq2) -
                       candidateWords[i].indexSeq1;

        if (diagonalVector[diagonal] == -1) {
          diagonalVector[diagonal] = candidateWords[i].indexSeq1;
          diagonalWordVector[diagonal] = candidateWords[i].wordSize;
        } else {
          // If the word is within distance of the previous same diagonal word
          if ((diagonalVector[diagonal] + distance >=
               candidateWords[i].indexSeq1)) {
            allowExtension = true;
          }

          // If the word is within distance of the next same diagonal word
          if ((i < candidateWords.size() - 1) && !allowExtension) {
            if (candidateWords[i].indexSeq1 + distance >=
                candidateWords[i + 1].indexSeq1) {
              allowExtension = true;
            }
          }

          diagonalVector[diagonal] = candidateWords[i].indexSeq1;
          diagonalWordVector[diagonal] = candidateWords[i].wordSize;
        }

        //<------------- Expand word ------------->
        ScoreSystemType score = candidateWords[i].score;
        bool expansion = false;
        if (allowExtension) {
          // Expand left
          int leftIndex = candidateWords[i].indexSeq1 - 1;
          if (leftIndex >= 0) { // We can expand left
            candidateWords[i].indexSeq1--;
            candidateWords[i].indexSeq2--;
            candidateWords[i].wordSize++;

            ScoreSystemType Similarity;

            if (candidateWords[i].indexSeq2 < 0) {
              Similarity = Mismatch;
            } else {
              Similarity = BaseType::match(Seq1[leftIndex],
                                           Seq2[candidateWords[i].indexSeq2])
                               ? Match
                               : Mismatch;
            }

            // A detremental expansion so revert
            if (candidateWords[i].score + Similarity <
                    candidateWords[i].score ||
                Similarity < 0) {
              candidateWords[i].indexSeq1++;
              candidateWords[i].indexSeq2++;
              candidateWords[i].wordSize--;
            } else {
              score += Similarity;
              expansion = true;
              anyExpansion = true;
            }
          }

          // Expand right
          int rightIndex =
              candidateWords[i].indexSeq1 +
              candidateWords[i].wordSize; // Plus wordSize as we index points
                                          // towards start of word
          if (rightIndex < Seq1.size()) {
            candidateWords[i].wordSize++;

            ScoreSystemType Similarity;

            if (candidateWords[i].indexSeq2 + candidateWords[i].wordSize >=
                Seq2.size()) {
              Similarity = Mismatch;
            } else {
              Similarity = BaseType::match(Seq1[rightIndex],
                                           Seq2[candidateWords[i].indexSeq2 +
                                                candidateWords[i].wordSize - 1])
                               ? Match
                               : Mismatch;
            }

            // A detremental expansion so revert
            if (candidateWords[i].score + Similarity <
                    candidateWords[i].score ||
                Similarity < 0) {
              candidateWords[i].wordSize--;
            } else {
              score += Similarity;
              expansion = true;
              anyExpansion = true;
            }
          }
        }

        // Store word that is complete (cannot be expanded further)
        // And is above a certain threshold
        // These words will be used for creating gapped alignments
        if ((!expansion) && (allowExtension) && (!candidateWords[i].saved) &&
            candidateWords[i].score > gapExpansionThreshold) {
          completeWords.push_back(candidateWords[i]);
          candidateWords[i].saved = true;
        }

        // Expanded seeds that improve upon alignment are saved for further
        // expansion
        if (score > candidateWords[i].score) //&& score > prevBestScore)
        {
          candidateWords[i].score = score;
          // prevBestScore = score;
        }
      }

      // No expansions recently so leave with our complete words
      if (!anyExpansion) {
        if (completeWords.size() == 0) {
          candidateWord optimal;
          int maxScore = -9999999;

          for (int i = 0; i < candidateWords.size(); i++) {

            if (candidateWords[i].score > maxScore) {
              maxScore = candidateWords[i].score;
              optimal = candidateWords[i];
            }
          }

          completeWords.clear();
          completeWords.push_back(optimal);
          break;
        }

        break;
      }
    }

    //<----- Longest Increasing Subsequence ----->
    // Find the longest increasing subsequence of the complete words
    // Important for smaller search space and overall score/usage of words

    std::sort(completeWords.begin(), completeWords.end());

    if (completeWords.size() > 1) {
      // Denoting X as completeWords
      int P[completeWords.size()];
      int M[completeWords.size() + 1];

      int l = 0;
      int lo = 0;
      int hi = 0;
      for (int i = 0; i < completeWords.size(); i++) {
        // Binary search for the largest positive j <= l
        // Such that X[M[j]] <= X[i]
        lo = 1;
        hi = l;
        while (lo <= hi) {
          int mid = ceil((lo + hi) / 2);
          if (completeWords[M[mid]] <= completeWords[i]) {
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
        newList[i] = completeWords[k];
        k = P[k];
      }

      completeWords.clear();
      completeWords = newList;
    }

    if (completeWords.size() == 0) {
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

    // Build alignment
    if (switched) {
      // ContainerType temp = Seq1;
      // Seq1 = Seq2;
      // Seq2 = temp;
      std::swap(Seq1, Seq2);

      for (int i = 0; i < completeWords.size(); i++) {
        int tempIndex = completeWords[i].indexSeq2;
        completeWords[i].indexSeq2 = completeWords[i].indexSeq1;
        completeWords[i].indexSeq1 = tempIndex;
      }
    }

    // For each complete HSP align with one another
    // For each word
    int i = 0;
    bool breakOut = false;
    for (i = 0; i < completeWords.size(); i++) {
      // Align word
      alignWord(completeWords[i], Seq1, Seq2, Result);

      // Bridge gap to word ignoring overlapping words
      int j = i + 1;
      while (true) {
        if (j < completeWords.size()) {
          // Overlapping word so look to bridge to next word if possible
          if (completeWords[i].indexSeq1 + completeWords[i].wordSize >
                  completeWords[j].indexSeq1 ||
              completeWords[i].indexSeq2 + completeWords[i].wordSize >
                  completeWords[j].indexSeq2) {
            j++;
            continue;
          }
          NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy>
              bridge(Scoring, BaseType::getMatchOperation());
          ArrayView<ContainerType> seq1Sub(Seq1);
          ArrayView<ContainerType> seq2Sub(Seq2);
          seq1Sub.sliceWindow(completeWords[i].indexSeq1 +
                                  completeWords[i].wordSize,
                              completeWords[j].indexSeq1);
          seq2Sub.sliceWindow(completeWords[i].indexSeq2 +
                                  completeWords[i].wordSize,
                              completeWords[j].indexSeq2);

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

    if (i == 0) {
      BaseType::forceGlobal(
          Seq1, Seq2, Result, completeWords[0].indexSeq1,
          completeWords[0].indexSeq2,
          completeWords[0].indexSeq1 + completeWords[0].wordSize,
          completeWords[0].indexSeq2 + completeWords[0].wordSize);
    } else {
      // Force global alignment
      BaseType::forceGlobal(
          Seq1, Seq2, Result, completeWords[0].indexSeq1,
          completeWords[0].indexSeq2,
          completeWords[i - 1].indexSeq1 + completeWords[i - 1].wordSize,
          completeWords[i - 1].indexSeq2 + completeWords[i - 1].wordSize);
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

  GappedBLATSA() : BaseType(getDefaultScoring(), nullptr) {}

  GappedBLATSA(ScoringSystem Scoring, MatchFnTy Match)
      : BaseType(Scoring, Match) {}

  virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1,
                                                  ContainerType &Seq2) {
    AlignedSequence<Ty, Blank> Result;
    buildAlignment(Seq1, Seq2, Result);
    return Result;
  }
};
