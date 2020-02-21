#include <vector>
#include <queue>
template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
class GappedBLATSA : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>
{
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
	const ScoreSystemType Mismatch = AllowMismatch
										 ? Scoring.getMismatchPenalty()
										 : std::numeric_limits<ScoreSystemType>::min();

	struct candidateWord
	{
		int score;	 //Score of the word
		int indexSeq1; //Index of the word on sequence 1
		int indexSeq2; //Index of the word on sequence 2
		int wordSize;
		int maxScore; //Maximum score achieved by the alignment
		bool saved = false;

		candidateWord *operator=(const candidateWord b)
		{
			this->score = b.score;
			this->indexSeq1 = b.indexSeq1;
			this->indexSeq2 = b.indexSeq2;
			this->wordSize = b.wordSize;
			this->maxScore = b.maxScore;
			this->saved = b.saved;
			return this;
		}

		bool operator<(const candidateWord &b) const
		{
			return indexSeq1 < b.indexSeq1;
		}

		bool operator<=(const candidateWord &b) const
		{
			return indexSeq1 <= b.indexSeq1;
		}
	};

	//Build the alignment using the BLAST hueristic algorithm
	void buildAlignment(ContainerType &Seq1, ContainerType &Seq2, AlignedSequence<Ty, Blank> &Result)
	{

		auto &Data = Result.Data;
		size_t SizeSeq1 = Seq1.size();
		size_t SizeSeq2 = Seq2.size();

		//Create diagonal vector that stores index of most recent hit on that diagonal
		std::vector<int> diagonalVector;
		std::vector<int> diagonalWordVector; //Stores word size of most recent diagonal so we can ignore overlapping words

		//Determine appropriate word size based on the size of the smallest sequence
		if (SizeSeq1 <= SizeSeq2)
		{
			if (SizeSeq1 < 500)
			{
				wordSize = (int)std::log2(SizeSeq1);
			}
			else
			{
				wordSize = 11;
			}
		}
		else
		{
			if (SizeSeq2 < 500)
			{
				wordSize = (int)std::log2(SizeSeq2);
			}
			else
			{
				wordSize = 11;
			}

			std::swap(SizeSeq1, SizeSeq2);
		}

		diagonalVector.resize(SizeSeq1 + SizeSeq2);
		diagonalWordVector.resize(SizeSeq1 + SizeSeq2);
		std::fill(diagonalVector.begin(), diagonalVector.end(), -1);
		std::fill(diagonalWordVector.begin(), diagonalWordVector.end(), -1);

		distance = 4 * wordSize;

		//Prefer to calculate overlapping k-mers of shorter sequence
		//Switch the variables as Seq1 is assumed to be the shorter sequence
		//Must remember that a switch occurred for alignment ordering purposes
		bool switched = false;
		if (SizeSeq1 > SizeSeq2)
		{
			ContainerType temp = Seq1;
			Seq1 = Seq2;
			Seq2 = temp;
			switched = true;
		}

		//Find words of sequence 1
		std::vector<candidateWord> words;
		std::vector<candidateWord> candidateWords;
		do
		{
			initialThreshold = Match * wordSize;

			for (int i = 0; i < Seq1.size() - wordSize + 1; i++)
			{
				candidateWord cWord;
				cWord.indexSeq1 = i;
				cWord.wordSize = wordSize;
				words.push_back(cWord);
			}

			//Compare k-mer words to non-overlapping k-mer words in Sequence 2 and a score for each comparison and delete alignments below a threshold
			for (int i = 0; i < words.size(); i++)
			{
				for (int j = 0; j < Seq2.size() - wordSize + 1; j = j + wordSize)
				{
					ScoreSystemType score = 0;

					for (int k = 0; k < wordSize; k++)
					{

						ScoreSystemType Similarity = BaseType::match(Seq1[words[i].indexSeq1 + k], Seq2[j + k]) ? Match : Mismatch;
						score += Similarity;
					}

					if (score >= initialThreshold)
					{
						candidateWord cWord;
						cWord = words[i];
						cWord.score = score;
						cWord.indexSeq2 = j;
						candidateWords.push_back(cWord);
					}
				}
			}

			if (candidateWords.size() == 0)
			{
				wordSize = (int)wordSize * 3 / 4;

				if (wordSize <= 0)
				{
					return;
				}
			}

			words.clear();
		} while (candidateWords.size() == 0);

		gapExpansionThreshold = (wordSize * Match) + (Match * log2(SizeSeq2)) - 5;
		//std::cout << "Gap expansion threshold: " << gapExpansionThreshold << std::endl;
		//std::cout << "Wordsize: " << wordSize << std::endl;
		//std::cout << "Distance: " << distance << std::endl;

		std::vector<candidateWord> completeWords;
		std::vector<candidateWord> newCandidates = candidateWords;

		//Begin expanding out the seeds in both directions
		while (candidateWords.size() >= 1)
		{
			std::vector<candidateWord> pushedWords;

			bool anyExpansion = false;
			for (int i = 0; i < candidateWords.size(); i++)
			{

				//Only expand if within distance A of the previous hit
				bool allowExtension = false;
				int diagonal = (SizeSeq1 + newCandidates[i].indexSeq2) - newCandidates[i].indexSeq1;

				if (diagonalVector[diagonal] == -1)
				{
					diagonalVector[diagonal] = newCandidates[i].indexSeq1;
					diagonalWordVector[diagonal] = newCandidates[i].wordSize;
				}
				else
				{
					//If the word is within distance of the previous same diagonal word
					if ((diagonalVector[diagonal] + distance >= newCandidates[i].indexSeq1))
					{
						allowExtension = true;
					}

					//If the word is within distance of the next same diagonal word
					if ((i < candidateWords.size() - 1) && !allowExtension)
					{
						if (newCandidates[i].indexSeq1 + distance >= newCandidates[i + 1].indexSeq1)
						{
							allowExtension = true;
						}
					}

					diagonalVector[diagonal] = newCandidates[i].indexSeq1;
					diagonalWordVector[diagonal] = newCandidates[i].wordSize;
				}

				//<------------- Expand word ------------->
				ScoreSystemType score = newCandidates[i].score;
				bool expansion = false;
				if (allowExtension)
				{
					//Expand left
					int leftIndex = newCandidates[i].indexSeq1 - 1;
					if (leftIndex >= 0)
					{ //We can expand left
						char leftChar = Seq1[leftIndex];
						newCandidates[i].indexSeq1--;
						newCandidates[i].indexSeq2--;
						newCandidates[i].wordSize++;

						ScoreSystemType Similarity;

						if (newCandidates[i].indexSeq2 < 0)
						{
							Similarity = Mismatch;
						}
						else
						{
							Similarity = BaseType::match(Seq1[leftIndex], Seq2[newCandidates[i].indexSeq2]) ? Match : Mismatch;
						}

						//A detremental expansion so revert
						if (newCandidates[i].score + Similarity < newCandidates[i].score)
						{
							newCandidates[i].indexSeq1++;
							newCandidates[i].indexSeq2++;
							newCandidates[i].wordSize--;
						}
						else
						{
							score += Similarity;
							expansion = true;
							anyExpansion = true;
						}
					}

					//Expand right
					int rightIndex = newCandidates[i].indexSeq1 + newCandidates[i].wordSize; //Plus wordSize as we index points towards start of word
					if (rightIndex <= Seq1.size())
					{
						char rightChar = Seq1[rightIndex];
						newCandidates[i].wordSize++;

						ScoreSystemType Similarity;

						if (newCandidates[i].indexSeq2 + newCandidates[i].wordSize > Seq2.size())
						{
							Similarity = Mismatch;
						}
						else
						{
							Similarity = BaseType::match(Seq1[rightIndex], Seq2[newCandidates[i].indexSeq2 + newCandidates[i].wordSize - 1]) ? Match : Mismatch;
						}

						//A detremental expansion so revert
						if (newCandidates[i].score + Similarity < newCandidates[i].score)
						{
							newCandidates[i].wordSize--;
						}
						else
						{
							score += Similarity;
							expansion = true;
							anyExpansion = true;
						}
					}
				}

				//Store word that is complete (cannot be expanded further)
				//And is above a certain threshold
				//These words will be used for creating gapped alignments
				if ((!expansion) && (allowExtension) && (!newCandidates[i].saved) && newCandidates[i].score > gapExpansionThreshold)
				{
					completeWords.push_back(newCandidates[i]);
					newCandidates[i].saved = true;
				}

				//Expanded seeds that improve upon alignment are saved for further expansion
				if (score > newCandidates[i].score) //&& score > prevBestScore)
				{
					newCandidates[i].score = score;
					//prevBestScore = score;
				}
			}

			//No expansions recently so leave with our complete words
			if (!anyExpansion)
			{
				if (completeWords.size() == 0)
				{
					candidateWord optimal;
					int maxScore = -9999999;

					for (int i = 0; i < candidateWords.size(); i++)
					{

						if (candidateWords[i].score > maxScore)
						{
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
		//Find the longest increasing subsequence of the complete words
		//Important for smaller search space and overall score/usage of words
		if (completeWords.size() != 1)
		{
			//Denoting X as completeWords
			int P[completeWords.size()];
			int M[completeWords.size() + 1];

			int l = 0;
			int lo = 0;
			int hi = 0;
			for (int i = 0; i < completeWords.size() - 1; i++)
			{
				//Binary search for the largest positive j <= l
				//Such that X[M[j]] <= X[i]
				lo = 1;
				hi = l;
				while (lo <= hi)
				{
					int mid = ceil((lo + hi) / 2);
					if (completeWords[M[mid]] <= completeWords[i])
					{
						lo = mid + 1;
					}
					else
					{
						hi = mid - 1;
					}
				}

				//After searching, lo is 1 greater than the length
				//of the longest prefix of X[i]
				int newL = lo;

				//The predecessor of X[i] is the last index of
				//the subsequence of length newL - 1
				P[i] = M[newL - 1];
				M[newL] = i;

				//If we found a subsequence longer than any we've foudn yet
				if (newL > l)
				{
					l = newL;
				}
			}

			//Reconstruct the longest increasing subsequence
			std::vector<candidateWord> newList;
			newList.resize(l);
			int k = M[l];
			for (int i = l - 1; i >= 0; i--)
			{
				newList[i] = completeWords[k];
				k = P[k];
			}

			completeWords.clear();
			completeWords = newList;
		}

		//Build alignment
		if (switched)
		{
			ContainerType temp = Seq1;
			Seq1 = Seq2;
			Seq2 = temp;

			for (int i = 0; i < completeWords.size(); i++)
			{
				int tempIndex = completeWords[i].indexSeq2;
				completeWords[i].indexSeq2 = completeWords[i].indexSeq1;
				completeWords[i].indexSeq1 = tempIndex;
			}
		}

		//For each complete HSP align with one another
		int recentIndexSeq1 = 0;
		int recentIndexSeq2 = 0;
		bool allowAlign = true;
		candidateWord recentWord;

		for (int i = 0; i < completeWords.size(); i++)
		{
			//Align HSP
			if (allowAlign)
			{
				for (int j = 0; j < completeWords[i].wordSize; j++)
				{
					bool isValidMatch = BaseType::match(Seq1[completeWords[i].indexSeq1 + j], Seq2[completeWords[i].indexSeq2 + j]);

					Data.push_back(
						typename BaseType::EntryType(Seq1[completeWords[i].indexSeq1 + j], Seq2[completeWords[i].indexSeq2 + j], isValidMatch));
				}

				recentWord = completeWords[i];
			}

			//Gap align towards next HSP
			if (i < completeWords.size() - 1)
			{
				if (completeWords[i + 1].indexSeq1 < recentWord.indexSeq1 + recentWord.wordSize || completeWords[i + 1].indexSeq2 < recentWord.indexSeq2 + recentWord.wordSize)
				{
					allowAlign = false;
					continue;
				}

				allowAlign = true;

				int indexSeq1 = recentWord.indexSeq1 + recentWord.wordSize;
				int indexSecondSeq1 = completeWords[i + 1].indexSeq1;
				int indexSeq2 = recentWord.indexSeq2 + recentWord.wordSize;
				int indexSecondSeq2 = completeWords[i + 1].indexSeq2;
				NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> SA(Scoring, BaseType::getMatchOperation());

				ArrayView<ContainerType> seq1Sub(Seq1);
				seq1Sub.sliceWindow(indexSeq1, indexSecondSeq1);
				ArrayView<ContainerType> seq2Sub(Seq2);
				seq2Sub.sliceWindow(indexSeq2, indexSecondSeq2);

				AlignedSequence<Ty, Blank> NWAlignment = SA.getAlignment(seq1Sub, seq2Sub);
				Data.insert(Data.end(), NWAlignment.begin(), NWAlignment.end());
			}
		}
	}

  public:
	static ScoringSystem
	getDefaultScoring()
	{
		return ScoringSystem(-1, 2, -1);
	}

	GappedBLATSA() : BaseType(getDefaultScoring(), nullptr) {}

	GappedBLATSA(ScoringSystem Scoring, MatchFnTy Match)
		: BaseType(Scoring, Match) {}

	virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1, ContainerType &Seq2)
	{
		AlignedSequence<Ty, Blank> Result;
		buildAlignment(Seq1, Seq2, Result);
		return Result;
	}
};
