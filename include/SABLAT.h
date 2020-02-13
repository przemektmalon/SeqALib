#include <vector>

template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
class BLATSA : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>
{
  private:
	int wordSize;
	int initialThreshold;

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
		int maxScore;		//Maximum score achieved by the alignment
		ContainerType word; //The word itself

		candidateWord *operator=(const candidateWord b)
		{
			this->score = b.score;
			this->indexSeq1 = b.indexSeq1;
			this->indexSeq2 = b.indexSeq2;
			this->wordSize = b.wordSize;
			this->maxScore = b.maxScore;
			this->word = b.word;
			return this;
		}
	};

	//Build the alignment using the BLAST hueristic algorithm
	void buildAlignment(ContainerType &Seq1, ContainerType &Seq2, AlignedSequence<Ty, Blank> &Result)
	{

		auto &Data = Result.Data;
		const size_t SizeSeq1 = Seq1.size();
		const size_t SizeSeq2 = Seq2.size();

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
		}

		//Prefer to calculate overlapping k-mers of shorter sequence
		//Switch the variables as Seq1 is assumed to be the shorter sequence
		//Must remember that a switch occurred for alignment ordering purposes
		bool switched = false;
		if (Seq1.size() > Seq2.size())
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
			//For near-perfect matches use Match*k-1
			initialThreshold = Match * (wordSize - 1);

			for (int i = 0; i < Seq1.size() - wordSize + 1; i++)
			{
				ContainerType word = Seq1.substr(i, wordSize);
				candidateWord cWord;
				cWord.word = word;
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

						char wordChar = words[i].word[k];
						ScoreSystemType Similarity = wordChar == Seq2[j + k] ? Match : Mismatch;

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

		//Begin expanding out the seeds in both directions
		while (candidateWords.size() >= 1)
		{
			std::vector<candidateWord> newCandidates = candidateWords;
			std::vector<candidateWord> pushedWords;
			for (int i = 0; i < newCandidates.size(); i++)
			{
				//<------------- Merge overlapping words ------------->
				if (i > 0)
				{
					candidateWord firstWord = newCandidates[i - 1];
					candidateWord secondWord = newCandidates[i];

					//The position at which the first word ends is greater than the beginning of the second word
					//This means they are overlapping so merge them together
					if ((firstWord.indexSeq1 + firstWord.wordSize > secondWord.indexSeq1) && (firstWord.indexSeq2 + firstWord.wordSize > secondWord.indexSeq2) && (firstWord.indexSeq1!=secondWord.indexSeq1))
					{

						int overlapSize = firstWord.indexSeq1 + firstWord.wordSize - secondWord.indexSeq1;
						int overlapSeq2Size = firstWord.indexSeq2 + firstWord.wordSize - secondWord.indexSeq2;
						if (overlapSize != overlapSeq2Size) {
							i++;
							continue;
						}

						candidateWord mergedWord = firstWord;
						mergedWord.word = mergedWord.word + secondWord.word.substr(overlapSize);

						//re-calculate score, however, we already know word 1 score so to reduce time we only calculate score of the word 2 that is being added to
						//word 1. Must be done since imperfect matches prevents us from simply multiplying a match score with the number of characters being added
						//to word 1
						int mergedWordSize = mergedWord.wordSize + secondWord.wordSize - (firstWord.indexSeq1 + firstWord.wordSize - secondWord.indexSeq1);
						for (int j = mergedWord.wordSize; j < mergedWordSize; j++)
						{
							char wordChar = mergedWord.word[j];
							ScoreSystemType Similarity = wordChar == Seq2[mergedWord.indexSeq2 + mergedWord.wordSize] ? Match : Mismatch;
							mergedWord.score += Similarity;
						}

						mergedWord.wordSize = mergedWordSize;
						newCandidates[i - 1] = mergedWord;
						newCandidates.erase(newCandidates.begin() + i);
						i--;
						continue;
					}
					//The two words are on the verge of joining so merging is simple addition
					//Separate else if for efficiency's sake
					else if (firstWord.indexSeq1 + firstWord.wordSize == secondWord.indexSeq1)
					{
						candidateWord mergedWord = firstWord;
						mergedWord.word = mergedWord.word + secondWord.word;
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
				//Expand left
				int leftIndex = newCandidates[i].indexSeq1 - 1;
				if (leftIndex >= 0)
				{ //We can expand left
					char leftChar = Seq1[leftIndex];
					newCandidates[i].word = leftChar + newCandidates[i].word;
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
						Similarity = leftChar == Seq2[newCandidates[i].indexSeq2] ? Match : Mismatch;
					}

					//A detremental expansion so revert
					if (newCandidates[i].score + Similarity < newCandidates[i].score)
					{
						newCandidates[i].word = newCandidates[i].word.substr(1, newCandidates[i].word.size());
						newCandidates[i].indexSeq1++;
						newCandidates[i].indexSeq2++;
						newCandidates[i].wordSize--;
					}
					else
					{
						score += Similarity;
					}
				}

				//Expand right
				int rightIndex = newCandidates[i].indexSeq1 + newCandidates[i].wordSize; //Plus wordSize as we index points towards start of word
				if (rightIndex <= Seq1.size())
				{
					char rightChar = Seq1[rightIndex];
					newCandidates[i].word = newCandidates[i].word + rightChar;
					newCandidates[i].wordSize++;

					ScoreSystemType Similarity;

					if (newCandidates[i].indexSeq2 + newCandidates[i].wordSize > Seq2.size())
					{
						Similarity = Mismatch;
					}
					else
					{
						Similarity = rightChar == Seq2[newCandidates[i].indexSeq2 + newCandidates[i].wordSize - 1] ? Match : Mismatch;
					}

					//A detremental expansion so revert
					if (newCandidates[i].score + Similarity < newCandidates[i].score)
					{
						newCandidates[i].word = newCandidates[i].word.substr(0, newCandidates[i].word.size() - 1);
						newCandidates[i].wordSize--;
					}
					else
					{
						score += Similarity;
					}
				}

				//Expanded seeds that improve upon alignment are saved for further expansion
				if (score > newCandidates[i].score) //&& score > prevBestScore)
				{
					newCandidates[i].score = score;
					pushedWords.push_back(newCandidates[i]);
					//prevBestScore = score;
				}
			}

			//No alignments acceptable after most recent expansion
			//So search previous list of candidates to find highest scoring one
			if (pushedWords.size() == 0)
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

				candidateWords.clear();
				candidateWords.push_back(optimal);
				break;
			}
			candidateWords.clear();
			candidateWords = pushedWords;
		}

		//Build alignment
		int j = 0;
		if (switched)
		{
			ContainerType temp = Seq1;
			Seq1 = Seq2;
			Seq2 = temp;
			int tempIndex = candidateWords[0].indexSeq2;
			candidateWords[0].indexSeq2 = candidateWords[0].indexSeq1;
			candidateWords[0].indexSeq1 = tempIndex;
		}

		for (int i = candidateWords[0].indexSeq2; i < candidateWords[0].indexSeq2 + candidateWords[0].wordSize; i++)
		{

			bool isValidMatch = candidateWords[0].word[j] == Seq2[i];

			Data.push_back(
				typename BaseType::EntryType(candidateWords[0].word[j], Seq2[i], isValidMatch));

			j++;
		}
	}

  public:
	static ScoringSystem getDefaultScoring()
	{
		return ScoringSystem(-1, 2, -1);
	}

	BLATSA() : BaseType(getDefaultScoring(), nullptr) {}

	BLATSA(ScoringSystem Scoring, MatchFnTy Match)
		: BaseType(Scoring, Match) {}

	virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1, ContainerType &Seq2)
	{
		AlignedSequence<Ty, Blank> Result;
		buildAlignment(Seq1, Seq2, Result);
		return Result;
	}
};
