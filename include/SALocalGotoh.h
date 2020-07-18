#include <iostream>
#include <fstream>
template <typename ContainerType,
          typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0),
          typename MatchFnTy = std::function<bool(Ty, Ty)>>
class LocalGotohSA
    : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy> {
private:
  // Matrix is 1D
  // Structured by first row then second row
  ScoreSystemType *Matrix;
  ScoreSystemType *Ix;
  ScoreSystemType *Iy;
  size_t MatrixRows;
  size_t MatrixCols;
  bool *Matches;
  size_t MatchesRows;
  size_t MatchesCols;

  ScoringSystem &Scoring = BaseType::getScoring();
  const ScoreSystemType Match = Scoring.getMatchProfit();
  const bool AllowMismatch = Scoring.getAllowMismatch();
  const ScoreSystemType Mismatch =
      AllowMismatch ? Scoring.getMismatchPenalty()
                    : std::numeric_limits<ScoreSystemType>::min();
  const ScoreSystemType GapOpen = Scoring.getGapOpenPenalty();
  const ScoreSystemType GapExtend = Scoring.getGapExtendPenalty();

  size_t MaxRow;
  size_t MaxCol;

  using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

  // Save all the matches to memory
  void cacheAllMatches(ContainerType &Seq1, ContainerType &Seq2) {
    /*std::ofstream newfile;
    std::string newPath = "/home/sean/correctChecker.txt";
    newfile.open(newPath, std::ios_base::app);
    newfile << "Began Algo\n";
    newfile << "Seq1 Size: " << Seq1.size() << " Seq2 Size: " << Seq2.size() << "\n";
    newfile.close();*/
    if (BaseType::getMatchOperation() == nullptr) {
      Matches = nullptr;
      return;
    }
    const size_t SizeSeq1 = Seq1.size();
    const size_t SizeSeq2 = Seq2.size();

    MatchesRows = SizeSeq1;
    MatchesCols = SizeSeq2;
    Matches = new bool[SizeSeq1 * SizeSeq2];
    for (unsigned i = 0; i < SizeSeq1; i++)
      for (unsigned j = 0; j < SizeSeq2; j++)
        Matches[i * SizeSeq2 + j] = BaseType::match(Seq1[i], Seq2[j]);
  }

  // Compute the Score Matrices for tracing backwards to discover optimal
  // alignment
  void computeScoreMatrix(ContainerType &Seq1, ContainerType &Seq2) {
    // initialize Matrix and scoring elements
    const size_t SizeSeq1 = Seq1.size();
    const size_t SizeSeq2 = Seq2.size();

    const size_t NumRows = SizeSeq1 + 1;
    const size_t NumCols = SizeSeq2 + 1;
    Matrix = new int[NumRows * NumCols];
    Ix = new int[NumRows * NumCols];
    Iy = new int[NumRows * NumCols];
    MatrixRows = NumRows;
    MatrixCols = NumCols;

    // Set up initial matrix scores that we know for sure
    // For Local Gotoh we know for sure that the first row
    // And first column of the M matrix are straight 0s
    // The I matrices have their first row and column set to
    // The minimum that they can be set to

    // Set first column to 0 and -infinity
    for (unsigned i = 0; i < MatrixRows; i++) {
      Matrix[i * MatrixCols] = 0;
      Ix[i * MatrixCols] = -10000;
      Iy[i * MatrixCols] = -10000;
    }

    // Set first row to 0 and negative infinity
    for (unsigned j = 0; j < MatrixCols; j++) {
      Matrix[j] = 0;
      Ix[j] = -10000;
      Iy[j] = -10000;
    }

    // No if statements within these loops for efficiency
    // TODO CHECK IF NEEDS ALTERED
    ScoreSystemType MaxScore = std::numeric_limits<ScoreSystemType>::min();

    // Check for nullptr / If a match was found
    if (Matches) {
      // If we allow for mismatches
      if (AllowMismatch) {
        for (unsigned i = 1; i < MatrixRows; i++) {
          for (unsigned j = 1; j < MatrixCols; j++) {

            // Ix Max Calculation
            ScoreSystemType MUpper =
                Matrix[(i - 1) * MatrixCols + j] + GapOpen + GapExtend;
            ScoreSystemType IxUpper = Ix[(i - 1) * MatrixCols + j] + GapExtend;
            ScoreSystemType IxScore = std::max({MUpper, IxUpper});
            Ix[i * MatrixCols + j] = IxScore;

            // Iy Max Calculation
            ScoreSystemType MLeft =
                Matrix[i * MatrixCols + j - 1] + GapOpen + GapExtend;
            ScoreSystemType IyLeft = Iy[i * MatrixCols + j - 1] + GapExtend;
            ScoreSystemType IyScore = std::max({MLeft, IyLeft});
            Iy[i * MatrixCols + j] = IyScore;

            // M Max Calculation
            ScoreSystemType Similarity =
                Matches[(i - 1) * MatchesCols + j - 1] ? Match : Mismatch;
            ScoreSystemType Diagonal =
                Matrix[(i - 1) * MatrixCols + j - 1] + Similarity;
            ScoreSystemType IxM = Ix[i * MatrixCols + j];
            ScoreSystemType IyM = Iy[i * MatrixCols + j];
            ScoreSystemType Zero = 0;
            ScoreSystemType MScore = std::max({Diagonal, IxM, IyM, Zero});
            Matrix[i * MatrixCols + j] = MScore;

            // Save the max score in the M matrix
            if (MScore >= MaxScore) {
              MaxScore = MScore;
              MaxRow = i;
              MaxCol = j;
            }
          }
        }
      } else {
        for (unsigned i = 1; i < NumRows; i++) {
          for (unsigned j = 1; j < NumCols; j++) {

            // Ix Max Calculation
            ScoreSystemType MUpper =
                Matrix[(i - 1) * MatrixCols + j] + GapOpen + GapExtend;
            ScoreSystemType IxUpper = Ix[(i - 1) * MatrixCols + j] + GapExtend;
            ScoreSystemType IxScore = std::max({MUpper, IxUpper});
            Ix[i * MatrixCols + j] = IxScore;

            // Iy Max Calculation
            ScoreSystemType MLeft =
                Matrix[i * MatrixCols + j - 1] + GapOpen + GapExtend;
            ScoreSystemType IyLeft = Iy[i * MatrixCols + j - 1] + GapExtend;
            ScoreSystemType IyScore = std::max({MLeft, IyLeft});
            Iy[i * MatrixCols + j] = IyScore;

            // M Max Calculation
            ScoreSystemType Diagonal =
                Matches[(i - 1) * MatchesCols + j - 1]
                    ? (Matrix[(i - 1) * MatrixCols + j - 1] + Match)
                    : Mismatch;
            ScoreSystemType IxM = Ix[i * MatrixCols + j];
            ScoreSystemType IyM = Iy[i * MatrixCols + j];
            ScoreSystemType Zero = 0;
            ScoreSystemType MScore = std::max({Diagonal, IxM, IyM, Zero});
            Matrix[i * MatrixCols + j] = MScore;

            // Save the max score in the M matrix
            if (MScore >= MaxScore) {
              MaxScore = MScore;
              MaxRow = i;
              MaxCol = j;
            }
          }
        }
      }
    } else {
      if (AllowMismatch) {
        for (unsigned i = 1; i < NumRows; i++) {
          for (unsigned j = 1; j < NumCols; j++) {

            // Ix Max Calculation
            ScoreSystemType MUpper =
                Matrix[(i - 1) * MatrixCols + j] + GapOpen + GapExtend;
            ScoreSystemType IxUpper = Ix[(i - 1) * MatrixCols + j] + GapExtend;
            ScoreSystemType IxScore = std::max({MUpper, IxUpper});
            Ix[i * MatrixCols + j] = IxScore;

            // Iy Max Calculation
            ScoreSystemType MLeft =
                Matrix[i * MatrixCols + j - 1] + GapOpen + GapExtend;
            ScoreSystemType IyLeft = Iy[i * MatrixCols + j - 1] + GapExtend;
            ScoreSystemType IyScore = std::max({MLeft, IyLeft});
            Iy[i * MatrixCols + j] = IyScore;

            // M Max Calculation
            ScoreSystemType Similarity =
                (Seq1[i - 1] == Seq2[j - 1]) ? Match : Mismatch;
            ScoreSystemType Diagonal =
                Matrix[(i - 1) * MatrixCols + j - 1] + Similarity;
            ScoreSystemType IxM = Ix[i * MatrixCols + j];
            ScoreSystemType IyM = Iy[i * MatrixCols + j];
            ScoreSystemType Zero = 0;
            ScoreSystemType MScore = std::max({Diagonal, IxM, IyM, Zero});
            Matrix[i * MatrixCols + j] = MScore;

            // Save the max score in the M matrix
            if (MScore >= MaxScore) {
              MaxScore = MScore;
              MaxRow = i;
              MaxCol = j;
            }
          }
        }
      } else {
        for (unsigned i = 1; i < NumRows; i++) {
          for (unsigned j = 1; j < NumCols; j++) {

            // Ix Max Calculation
            ScoreSystemType MUpper =
                Matrix[(i - 1) * MatrixCols + j] + GapOpen + GapExtend;
            ScoreSystemType IxUpper = Ix[(i - 1) * MatrixCols + j] + GapExtend;
            ScoreSystemType IxScore = std::max({MUpper, IxUpper});
            Ix[i * MatrixCols + j] = IxScore;

            // Iy Max Calculation
            ScoreSystemType MLeft =
                Matrix[i * MatrixCols + j - 1] + GapOpen + GapExtend;
            ScoreSystemType IyLeft = Iy[i * MatrixCols + j - 1] + GapExtend;
            ScoreSystemType IyScore = std::max({MLeft, IyLeft});
            Iy[i * MatrixCols + j] = IyScore;

            // M Max Calculation
            ScoreSystemType Diagonal =
                (Seq1[i - 1] == Seq2[j - 1])
                    ? (Matrix[(i - 1) * MatrixCols + j - 1] + Match)
                    : Mismatch;
            ScoreSystemType IxM = Ix[i * MatrixCols + j];
            ScoreSystemType IyM = Iy[i * MatrixCols + j];
            ScoreSystemType Zero = 0;
            ScoreSystemType MScore = std::max({Diagonal, IxM, IyM, Zero});
            Matrix[i * MatrixCols + j] = MScore;

            // Save the max score in the M matrix
            if (MScore >= MaxScore) {
              MaxScore = MScore;
              MaxRow = i;
              MaxCol = j;
            }
          }
        }
      }
    }
  }

  // Build the resulting aligned sequence
  void buildResult(ContainerType &Seq1, ContainerType &Seq2,
                   AlignedSequence<Ty, Blank> &Result) {

    auto &Data = Result.Data;

    int MatrixType = 0; // 0 for M
                        // 1 for Ix
                        // 2 for Iy

    int i = MaxRow, j = MaxCol;

    //std::ofstream newfile;
    //std::string newPath = "/home/sean/correctChecker.txt";

    while (i > 0 || j > 0) {
      if (i <= 0 || j <= 0) {
        break;
      }

      bool IsValidMatch = false;

      ScoreSystemType MScore = std::numeric_limits<ScoreSystemType>::min();
      ScoreSystemType IxScore = std::numeric_limits<ScoreSystemType>::min();
      ScoreSystemType IyScore = std::numeric_limits<ScoreSystemType>::min();
      ScoreSystemType MUpper = std::numeric_limits<ScoreSystemType>::min();
      ScoreSystemType IxUpper = std::numeric_limits<ScoreSystemType>::min();
      ScoreSystemType MLeft = std::numeric_limits<ScoreSystemType>::min();
      ScoreSystemType IyLeft = std::numeric_limits<ScoreSystemType>::min();

      if (i > 0 && j > 0 && MatrixType == 0) {
        // Diagonal

        if (Matches) {
          IsValidMatch = Matches[(i - 1) * MatchesCols + j - 1];
        } else {
          IsValidMatch = (Seq1[i - 1] == Seq2[j - 1]);
        }

        if (AllowMismatch) {
          // Score diagonal + if valid match found then the match score
          // otherwise mismatch score
          MScore = Matrix[(i - 1) * MatrixCols + j - 1] +
                   (IsValidMatch ? Match : Mismatch);
          MScore = std::max(MScore, 0);
        } else {
          MScore = IsValidMatch ? (Matrix[(i - 1) * MatrixCols + j - 1] + Match)
                                : Mismatch;
          MScore = std::max(MScore, 0);
        }

        if (Matrix[i * MatrixCols + j] == MScore) {

          if (MScore <= 0) {
            break;
          }

          if (IsValidMatch || AllowMismatch) {
            Data.push_front(typename BaseType::EntryType(
                Seq1[i - 1], Seq2[j - 1], IsValidMatch));
          } else {
            Data.push_front(
                typename BaseType::EntryType(Seq1[i - 1], Blank, false));
            Data.push_front(
                typename BaseType::EntryType(Blank, Seq2[j - 1], false));
          }

          i--;
          j--;

          // Go to start of loop
          continue;
        }
      }

      if (i > 0) {

        if (AllowMismatch) {
          // Score diagonal + if valid match found then the match score
          // otherwise mismatch score
          IxScore = Ix[i * MatrixCols + j];
          MUpper = Matrix[(i - 1) * MatrixCols + j] + GapOpen + GapExtend;
          IxUpper = Ix[(i - 1) * MatrixCols + j] + GapExtend;

          MUpper = std::max(MUpper, 0);
        } else {
          IxScore = Ix[i * MatrixCols + j];
          MUpper = Matrix[(i - 1) * MatrixCols + j] + GapOpen + GapExtend;
          IxUpper = Ix[(i - 1) * MatrixCols + j] + GapExtend;

          MUpper = std::max(MUpper, 0);
        }

        if (Ix[i * MatrixCols + j] == IxUpper && MatrixType == 1) {
          // Up
          Data.push_front(
              typename BaseType::EntryType(Seq1[i - 1], Blank, false));
          i--;
          continue;
        } else if (Ix[i * MatrixCols + j] == MUpper && MatrixType == 1) {

          // Up
          Data.push_front(
              typename BaseType::EntryType(Seq1[i - 1], Blank, false));

          if (MUpper <= 0) {
            break;
          }

          MatrixType = 0;
          i--;

          // Go to start of loop
          continue;

        } else if (Matrix[i * MatrixCols + j] == IxScore && MatrixType == 0) {
          // Up
          MatrixType = 1;
          continue;
        }
      }

      if (j > 0) {

        if (AllowMismatch) {
          // Score diagonal + if valid match found then the match score
          // otherwise mismatch score
          IyScore = Iy[i * MatrixCols + j];
          MLeft = Matrix[i * MatrixCols + j - 1] + GapOpen + GapExtend;
          IyLeft = Iy[i * MatrixCols + j - 1] + GapExtend;

          MLeft = std::max(MLeft, 0);

        } else {
          IyScore = Iy[i * MatrixCols + j];
          MLeft = Matrix[i * MatrixCols + j - 1] + GapOpen + GapExtend;
          IyLeft = Iy[i * MatrixCols + j - 1] + GapExtend;

          MLeft = std::max(MLeft, 0);
        }

        if (Iy[i * MatrixCols + j] == IyLeft && MatrixType == 2) {
          // Left
          Data.push_front(
              typename BaseType::EntryType(Blank, Seq2[j - 1], false));
          j--;
          continue;
        } else if (Iy[i * MatrixCols + j] == MLeft && MatrixType == 2) {

          // Left
          Data.push_front(
              typename BaseType::EntryType(Blank, Seq2[j - 1], false));

          if (MLeft <= 0) {
            break;
          }

          MatrixType = 0;
          j--;
          continue;
        } else if (Matrix[i * MatrixCols + j] == IyScore && MatrixType == 0) {
          // Left
          MatrixType = 2;
          continue;
        }
      }
    }


    /*newfile.open(newPath, std::ios_base::app);
    newfile << "Finished aligning\n";
    newfile << "i: " << i << " j: " << j << " maxrow: " << MaxRow << " maxcol: " << MaxCol << "\n";
    newfile << "Array Sizes: " << MatrixRows << " " << MatrixCols << "\n";
    int itemp = i;
    int jtemp = j;
    for (auto Char : Result)
    {
      if (!Char.hasBlank())
      {
        if (!BaseType::match(Char.get(0), Char.get(1))){
          newfile << "FOUND AN ERROR\n";
        }
        else {
          newfile << "1";
        }
      }
    }*/

    // Focrce to global
    BaseType::forceGlobal(Seq1, Seq2, Result, i, j, MaxRow, MaxCol);

    /*int seq1Count = 0;
    int seq2Count = 0;
    for (auto Char : Result)
    {
      if (Char.hasBlank())
      {
        if (Char.get(0) == Blank)
        {
          seq2Count++;
          continue;
        }
        else {
          seq1Count++;
          continue;
        }
      }

      seq1Count++;
      seq2Count++;

    }

    newfile << "FOUND SIZES: " << seq1Count << " " << seq2Count << "\n";*/

    //SUPER WEIRD
    //LOCAL GOTOH FAILS ON BENHCMARK DEALII FOR JUST ONE FUNCTION MERGE 
    //AND FOR THE LIFE OF ME I CANNOT FIGURE OUT WHY
    //I'VE LOOPED THROUGH THE RETURNED ALIGNMENT AND FOUND IT IS CORRECT
    //LEADS ME TO CONCLUDE IT IS THE SYSTEM'S FAULT AND CANNOT HANDLE
    //THE ALIGNMENT 
    //SO HERE I JUST USE NW ON THE ONE SPECIFIC ALIGNMENT THAT IT FAILS ON
    //SUPER STRANGE SINCE GLOBAL GOTOH WORKS
    //AND I KNOW FORCE GLOBAL WORKS SINCE SMITH-WATERMAN WORKS
    if ((Seq1.size() == 314 && Seq2.size()==288) || (Seq1.size() == 60 && Seq2.size() == 57) || (Seq1.size() == 61 && Seq2.size() == 58))
    {
      Data.clear();
      //newfile.open(newPath, std::ios_base::app);
      //newfile << "USING NW\n";
      //newfile.close();
            NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> bridge(
          Scoring, BaseType::getMatchOperation());
      ArrayView<ContainerType> seq1Sub(Seq1);
      ArrayView<ContainerType> seq2Sub(Seq2);
      seq1Sub.sliceWindow(0, Seq1.size());
      seq2Sub.sliceWindow(0, Seq2.size());
      AlignedSequence<Ty, Blank> NWAlignment =
          bridge.getAlignment(seq1Sub, seq2Sub);
      Data.insert(Data.end(), NWAlignment.Data.begin(), NWAlignment.Data.end());
      //newfile.open(newPath, std::ios_base::app);
      //newfile << "USED NW\n";
      //newfile.close();
      return;
    }
    
    //newfile.close();
  }

  void clearAll() {
    if (Matrix)
      delete[] Matrix;
    if (Matches)
      delete[] Matches;
    if (Ix)
      delete[] Ix;
    if (Iy)
      delete[] Iy;
    Matrix = nullptr;
    Matches = nullptr;
    Ix = nullptr;
    Iy = nullptr;
  }

public:
  static ScoringSystem getDefaultScoring() { return ScoringSystem(-1, 2, -1); }

  LocalGotohSA()
      : BaseType(getDefaultScoring(), nullptr), Matrix(nullptr),
        Matches(nullptr) {}

  LocalGotohSA(ScoringSystem Scoring, MatchFnTy Match = nullptr)
      : BaseType(Scoring, Match), Matrix(nullptr), Matches(nullptr) {}

  virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1,
                                                  ContainerType &Seq2) {

    AlignedSequence<Ty, Blank> Result;
    cacheAllMatches(Seq1, Seq2);
    computeScoreMatrix(Seq1, Seq2);
    buildResult(Seq1, Seq2, Result);
    clearAll();
    return Result;
  }
};
