#include <queue>
#include <functional>

int tops = 0;
int score = 0;
int expanded = 0;
template <typename ContainerType, typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0), typename MatchFnTy = std::function<bool(Ty, Ty)>>
class FOGSAASA : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>
{
private:
    bool *Matches;
    size_t MatchesRows;
    size_t MatchesCols;
    int M;
    int N;
    size_t SizeSeq1;
    size_t SizeSeq2;

    ScoringSystem &Scoring = BaseType::getScoring();
    const ScoreSystemType Match = Scoring.getMatchProfit();
    const bool AllowMismatch = Scoring.getAllowMismatch();
    const ScoreSystemType Mismatch = AllowMismatch ? Scoring.getMismatchPenalty() : std::numeric_limits<ScoreSystemType>::min();
    const ScoreSystemType GapOpen = Scoring.getGapOpenPenalty();
    const ScoreSystemType GapExtend = Scoring.getGapExtendPenalty();

    using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

    enum AlignType
    {
        MM, // Match or Mismatch
        _G, // Gap in first sequence
        G_  // Gap in second sequence
    };

    class Node
    {
    public:
        Node() {}

        int P1 = 0;
        int P2 = 0;
        int presentScore = std::numeric_limits<int>::min();
        int Tmin = std::numeric_limits<int>::min();
        int Tmax = std::numeric_limits<int>::min();
        AlignType type;

        Node *operator=(const Node &b)
        {
            this->Tmax = b.Tmax;
            this->Tmin = b.Tmin;
            this->P1 = b.P1;
            this->P2 = b.P2;
            this->presentScore = b.presentScore;
            this->type = b.type;
            return this;
        }

        bool operator<(Node &b) const { return Tmax > b.Tmax; }

        // Hash sorts by Tmax, but on case where Tmax are the same
        // use Tmin
        bool operator>(const Node &b) const { return Tmin < b.Tmin; }

        void calculateScores(int newP1, int newP2, int M, int N, AlignType aType, ScoringSystem &scoreSystem)
        {

            const ScoreSystemType Gap = scoreSystem.getGapPenalty();
            const ScoreSystemType Match = scoreSystem.getMatchProfit();
            const bool AllowMismatch = scoreSystem.getAllowMismatch();
            const ScoreSystemType Mismatch = AllowMismatch ? scoreSystem.getMismatchPenalty() : std::numeric_limits<ScoreSystemType>::min();

            type = aType;
            P1 = newP1;
            P2 = newP2;

            // Calculate Future scores (Fmin, Fmax)
            int Fmax = 0;
            int Fmin = 0;
            int x1 = (N - P2);
            int x2 = (M - P1);
            if (x2 < x1)
            {
                Fmin = (x2 * Mismatch) + (Gap * (x1 - x2));
                Fmax = (x2 * Match) + (Gap * (x1 - x2));
            }
            else
            {
                Fmin = (x1 * Mismatch) + (Gap * (x2 - x1));
                Fmax = (x1 * Match) + (Gap * (x2 - x1));
            }

            int Ftemp = (x2 * Gap) + (x1 * Gap);

            if (!AllowMismatch)
            {
                Fmin = Ftemp;
            }

            // Add Present score to the Future score
            // obtaining the fitness score
            Tmin = presentScore + Fmin;
            Tmax = presentScore + Fmax;
        }
    };

    class HashPriorityQueue
    {
        int numBuckets;
        int min;
        int max;
        bool AllowMismatch = true;

        std::vector<std::priority_queue<Node, std::vector<Node>, std::greater<Node>>> *table;

    public:
        HashPriorityQueue(int min, int max, ScoringSystem &scoreSystem)
        {
            this->numBuckets = max - min;
            this->min = min;
            this->max = max;
            table = new std::vector<std::priority_queue<Node, std::vector<Node>, std::greater<Node>>>;
            table->resize(numBuckets);
            AllowMismatch = scoreSystem.getAllowMismatch();
        }

        int insertItem(int key, Node &node, int maxPointer)
        {
            int index = hashFunction(key);
            if (!checkIndex(index))
            {
                if (AllowMismatch)
                {
                    int a = key - (table->size() - index) - 1;
                    return a;
                }
                // This used to prevent bad indices
                // However when using no mismatches this should
                // Instead return the current maxPointer
                return maxPointer;
            }

            if (index >= table->size())
            {
                table->resize(numBuckets + 1 + (index - numBuckets));
                numBuckets = numBuckets + 1 + (index - numBuckets);
            }
            table->at(index).push(node);

            // First entry so maxPointer is set to the one being inserted
            if (maxPointer == min - 1)
            {
                return node.Tmax;
            }

            // If the one we are inserting is better than what we are currently
            // Pointing towards then return
            // A new pointer to that position in the hash
            if (node.Tmax > maxPointer)
            {
                maxPointer = node.Tmax;
            }

            return maxPointer;
        }

        // Get top from the hashed queue
        Node getTop(int key)
        {
            int index = hashFunction(key);

            tops++;

            if (!checkIndex(index))
            {
                Node a;
                return a;
            }

            if (table->at(index).size() == 0)
            {
                Node a;
                return a;
            }

            return table->at(index).top();
        }

        // Delete top item and return a new maxPointer
        int deleteTop(int key)
        {
            int index = hashFunction(key);

            if (!checkIndex(index))
            {
                int a = key - (table->size() - index) - 1;
                return a;
            }

            if (!table->at(index).empty())
            {
                table->at(index).pop();
            }

            int maxPointer = key;

            int i = index;
            while (i < table->size())
            {

                if (!table->at(i).empty())
                {
                    break;
                }

                i++;
                maxPointer--;
            }

            return maxPointer;
        }

        bool checkIndex(int index)
        {
            // Bad index so return false so maxPointer can be set to minimum
            // This should be caught and the current best path is returned
            if (index >= table->size() || index < 0)
            {
                return false;
            }

            return true;
        }

        int hashFunction(int key) { return max - key; }

        void deleteTable() { delete table; }
    };

    Node **c;
    HashPriorityQueue *hpqueue = nullptr;

    // Save all the matches to memory
    void cacheAllMatches(ContainerType &Seq1, ContainerType &Seq2)
    {
        if (BaseType::getMatchOperation() == nullptr)
        {
            Matches = nullptr;
            return;
        }
        SizeSeq1 = Seq1.size();
        SizeSeq2 = Seq2.size();

        MatchesRows, M = SizeSeq1;
        MatchesCols, N = SizeSeq2;
        Matches = new bool[SizeSeq1 * SizeSeq2];
        c = new Node *[(SizeSeq1 + 1) * (SizeSeq2 + 1) + 2];
        for (unsigned i = 0; i < SizeSeq1; i++)
            for (unsigned j = 0; j < SizeSeq2; j++)
                Matches[i * SizeSeq2 + j] = BaseType::match(Seq1[i], Seq2[j]);

        for (int i = 0; i < (SizeSeq1 + 1) * (SizeSeq2 + 1) + 2; i++)
        {
            c[i] = nullptr;
        }
    }

    // Build the resulting aligned sequence
    void buildAlignment(ContainerType &Seq1, ContainerType &Seq2,
                        AlignedSequence<Ty, Blank> &Result)
    {

        ScoringSystem &Scoring = BaseType::getScoring();
        const ScoreSystemType Gap = Scoring.getGapPenalty();
        const ScoreSystemType Match = Scoring.getMatchProfit();
        const bool AllowMismatch = Scoring.getAllowMismatch();
        const ScoreSystemType Mismatch = AllowMismatch ? Scoring.getMismatchPenalty() : std::numeric_limits<ScoreSystemType>::min();

        int P1 = 0;
        int P2 = 0;
        int optimal = std::numeric_limits<int>::min();

        int newTmax = 0;
        int maxTmax = 0;
        Node currentNode;
        currentNode.presentScore = 0;

        int Fmax = 0;
        int Fmin = 0;
        int x1 = (N - P2);
        int x2 = (M - P1);
        if (x2 < x1)
        {
            Fmin = (x2 * Mismatch) + (Gap * (x1 - x2));
            Fmax = (x2 * Match) + (Gap * (x1 - x2));
        }
        else
        {
            Fmin = (x1 * Mismatch) + (Gap * (x2 - x1));
            Fmax = (x1 * Match) + (Gap * (x2 - x1));
        }
        int Ftemp = (M * Gap) + (N * Gap);

        if (Fmin > Ftemp || !AllowMismatch)
        {
            Fmin = Ftemp;
        }

        currentNode.Tmax = Fmax;
        currentNode.Tmin = Fmin;

        int ml, sl, th;
        if (M > N)
        {
            ml = M;
            sl = N;
        }
        else
        {
            ml = N;
            sl = M;
        }
        th = ml * 30 / 100;

        int threshold = th * Match + (sl - th) * Mismatch + Gap * (ml - sl);

        int inserted = 0;

        hpqueue = new HashPriorityQueue(Fmin, Fmax, Scoring);
        int maxPointer = Fmin - 1; // Points towards top of hashed priority queue
                                   // Initially points towards nothing

        int pathend = 0;

        bool oneLoop = false;

        if (M != 0 && N != 0)
        {
            do
            {
                while (P1 <= (M - 1) || P2 <= (N - 1))
                {
                    auto t1 = std::chrono::high_resolution_clock::now();

                    if (!c[P1 * N + P2])
                    {
                        c[P1 * N + P2] = new Node();
                    }

                    if (currentNode.Tmax < c[P1 * N + P2]->Tmax)
                    {
                        // This can occur when the top of the queue has actually been done
                        // better beforehand this check is needed so an umpromising node
                        // popped from the queue isnt expanded upon Prune the current branch
                        currentNode = hpqueue->getTop(maxPointer);
                        maxPointer = hpqueue->deleteTop(maxPointer);

                        if ((maxPointer <= optimal || maxPointer == Fmin - 1) &&
                            oneLoop == true)
                        {
                            align(Seq1.size(), Seq2.size(), Seq1, Seq2, Result);
                            return;
                        }

                        P1 = currentNode.P1;
                        P2 = currentNode.P2;

                        continue;
                    }

                    if (currentNode.presentScore > c[P1 * N + P2]->presentScore)
                    {
                        *c[P1 * N + P2] = currentNode;
                    }

                    // Select the best child from the remaining children according to the
                    // Tmax Develop sequence by adding either Match, Mismatch, or Gap From
                    // each developed sequence, calculate fitness scores
                    bool IsValidMatch = false;
                    if (Matches)
                    {
                        IsValidMatch = Matches[P1 * N + P2];
                    }
                    else
                    {
                        IsValidMatch = (BaseType::match(Seq1[P1], Seq2[P2]));
                    }

                    Node MMNode;
                    Node _GNode;
                    Node G_Node;
                    Node childNode;

                    ScoreSystemType Similarity;
                    if (AllowMismatch)
                    {
                        Similarity = IsValidMatch ? Match : Mismatch;
                        MMNode.presentScore = currentNode.presentScore + Similarity;

                        // Integer Underflow
                        if (currentNode.presentScore < 0 && Similarity < 0 &&
                            MMNode.presentScore > 0)
                        {
                            MMNode.presentScore = std::numeric_limits<int>::min() + 1;
                        }
                    }
                    else
                    {
                        if (IsValidMatch)
                        {
                            Similarity = Match;
                            MMNode.presentScore = currentNode.presentScore + Similarity;
                            // Integer Underflow
                            if (currentNode.presentScore < 0 && Similarity < 0 &&
                                MMNode.presentScore > 0)
                            {
                                MMNode.presentScore = std::numeric_limits<int>::min() + 1;
                            }
                        }
                        else
                        {
                            Similarity = Mismatch;
                            MMNode.presentScore = Similarity;
                        }
                    }

                    _GNode.presentScore = currentNode.presentScore + Gap;
                    G_Node.presentScore = currentNode.presentScore + Gap;

                    // Integer Underflow
                    if (currentNode.presentScore < 0 && Gap < 0 &&
                        _GNode.presentScore > 0)
                    {
                        _GNode.presentScore = std::numeric_limits<int>::min() + 1;
                        G_Node.presentScore = std::numeric_limits<int>::min() + 1;
                    }

                    MMNode.calculateScores(P1 + 1, P2 + 1, M, N, AlignType::MM, Scoring);
                    _GNode.calculateScores(P1, P2 + 1, M, N, AlignType::_G, Scoring);
                    G_Node.calculateScores(P1 + 1, P2, M, N, AlignType::G_, Scoring);

                    // std::cout << "(" << P1 << "," << P2 << ")" << std::endl;
                    // std::cout << "[" << currentNode.Tmax << "," << currentNode.Tmin <<
                    // "]" << std::endl;
                    if (P1 > (M - 1))
                    {
                        MMNode.Tmax = std::numeric_limits<int>::min();
                        G_Node.Tmax = std::numeric_limits<int>::min();
                    }
                    if (P2 > (N - 1))
                    {
                        MMNode.Tmax = std::numeric_limits<int>::min();
                        _GNode.Tmax = std::numeric_limits<int>::min();
                    }

                    if (MMNode.Tmax >= std::max(_GNode.Tmax, G_Node.Tmax))
                    {
                        childNode = MMNode;
                        P1++;
                        P2++;
                        inserted += 2;
                        // maxPointer = hpqueue->insertItem(G_Node.Tmax, G_Node,
                        // maxPointer); maxPointer = hpqueue->insertItem(_GNode.Tmax,
                        // _GNode, maxPointer);

                        if (_GNode.Tmax >= Fmin && _GNode.Tmax > G_Node.Tmax)
                        {
                            inserted++;
                            maxPointer = hpqueue->insertItem(_GNode.Tmax, _GNode, maxPointer);
                            // maxPointer = hpqueue->insertItem(G_Node.Tmax, G_Node,
                            // maxPointer);
                        }
                        else if (G_Node.Tmax >= Fmin && G_Node.Tmax >= _GNode.Tmax)
                        {
                            inserted++;
                            maxPointer = hpqueue->insertItem(G_Node.Tmax, G_Node, maxPointer);
                            // maxPointer = hpqueue->insertItem(_GNode.Tmax, _GNode,
                            // maxPointer);
                        }
                    }
                    else if (G_Node.Tmax > std::max(MMNode.Tmax, _GNode.Tmax))
                    {
                        childNode = G_Node;
                        P1++;

                        // maxPointer = hpqueue->insertItem(_GNode.Tmax, _GNode,
                        // maxPointer); maxPointer = hpqueue->insertItem(MMNode.Tmax,
                        // MMNode, maxPointer);

                        if (_GNode.Tmax >= Fmin && _GNode.Tmax > MMNode.Tmax)
                        {
                            inserted++;
                            maxPointer = hpqueue->insertItem(_GNode.Tmax, _GNode, maxPointer);
                            // maxPointer = hpqueue->insertItem(MMNode.Tmax, MMNode,
                            // maxPointer);
                        }
                        else if (MMNode.Tmax >= Fmin && MMNode.Tmax >= _GNode.Tmax)
                        {
                            inserted++;
                            maxPointer = hpqueue->insertItem(MMNode.Tmax, MMNode, maxPointer);
                            // maxPointer = hpqueue->insertItem(_GNode.Tmax, _GNode,
                            // maxPointer);
                        }
                    }
                    else
                    {
                        childNode = _GNode;
                        P2++;

                        // maxPointer = hpqueue->insertItem(G_Node.Tmax, G_Node,
                        // maxPointer); maxPointer = hpqueue->insertItem(MMNode.Tmax,
                        // MMNode, maxPointer);

                        if (G_Node.Tmax >= Fmin && G_Node.Tmax > MMNode.Tmax)
                        {
                            inserted++;
                            maxPointer = hpqueue->insertItem(G_Node.Tmax, G_Node, maxPointer);
                            // maxPointer = hpqueue->insertItem(MMNode.Tmax, MMNode,
                            // maxPointer);
                        }
                        else if (MMNode.Tmax >= Fmin && MMNode.Tmax >= G_Node.Tmax)
                        {
                            inserted++;
                            // maxPointer = hpqueue->insertItem(G_Node.Tmax, G_Node,
                            // maxPointer);
                            maxPointer = hpqueue->insertItem(MMNode.Tmax, MMNode, maxPointer);
                        }
                    }
                    expanded++;

                    if (!c[P1 * N + P2])
                    {
                        c[P1 * N + P2] = new Node;
                    }

                    if (childNode.presentScore <= c[P1 * N + P2]->presentScore)
                    {

                        // Prune the current branch
                        childNode = hpqueue->getTop(maxPointer);
                        maxPointer = hpqueue->deleteTop(maxPointer);

                        if (maxPointer <= optimal || maxPointer == Fmin - 1 && oneLoop)
                        {
                            align(Seq1.size(), Seq2.size(), Seq1, Seq2, Result);
                            return;
                        }

                        P1 = childNode.P1;
                        P2 = childNode.P2;
                    }
                    else
                    {
                        *c[P1 * N + P2] = childNode;
                        if (childNode.Tmax < optimal)
                        {
                            // Prune the current branch
                            childNode = hpqueue->getTop(maxPointer);
                            maxPointer = hpqueue->deleteTop(maxPointer);

                            if (maxPointer <= optimal || maxPointer == Fmin - 1 && oneLoop)
                            {
                                align(Seq1.size(), Seq2.size(), Seq1, Seq2, Result);
                                return;
                            }
                            P1 = childNode.P1;
                            P2 = childNode.P2;
                        }
                    }

                    currentNode = childNode;

                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto duration =
                        std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                            .count();
                    // std::cout << "One loop: " << duration << std::endl;
                }

                oneLoop = true;

                if (!c[P1 * N + P2])
                {
                    c[P1 * N + P2] = new Node;
                    *c[P1 * N + P2] = currentNode;
                }

                if (c[P1 * N + P2]->presentScore > optimal)
                {
                    optimal = c[P1 * N + P2]->presentScore;
                }

                // If the greedy search has achieved the highest possible score
                // then align
                if (optimal == Fmax)
                {
                    align(Seq1.size(), Seq2.size(), Seq1, Seq2, Result);
                    return;
                }

                currentNode = hpqueue->getTop(maxPointer);
                maxPointer = hpqueue->deleteTop(maxPointer);

                if (maxPointer <= optimal || maxPointer == Fmin - 1)
                {

                    align(Seq1.size(), Seq2.size(), Seq1, Seq2, Result);
                    return;
                }

                P1 = currentNode.P1;
                P2 = currentNode.P2;
                newTmax = currentNode.Tmax;

                // New similarity calculations
                int cpl;
                if (P1 > P2)
                {
                    cpl = P1;
                }
                else
                {
                    cpl = P2;
                }

                // If top most node has Tmax so less than 30% similarity then
                // end the process and report approximate score
                if (((cpl > (70 * ml / 100)) && (optimal < threshold)) ||
                    (newTmax < threshold))
                {
                    //std::cout << "Exiting because lack of similarity" << std::endl;
                    align(Seq1.size(), Seq2.size(), Seq1, Seq2, Result);
                    //std::cout << "Expanded nodes: " << expanded << std::endl;
                    return;
                }

            } while (optimal < newTmax);
        }

        align(Seq1.size(), Seq2.size(), Seq1, Seq2, Result);
    }

    void align(int P1, int P2, ContainerType &Seq1, ContainerType &Seq2,
               AlignedSequence<Ty, Blank> &Result)
    {
        auto &Data = Result.Data;

        //<----- Backtrace from P1 and P2 using c[] align types ----->

        score = c[P1 * N + P2]->presentScore;

        while (P1 > 0 || P2 > 0)
        {

            AlignType type = c[P1 * N + P2]->type;

            // Match or mismatch
            if (type == AlignType::MM)
            {

                bool IsValidMatch = false;
                if (Matches)
                {
                    IsValidMatch = Matches[(P1 - 1) * N + (P2 - 1)];
                }
                else
                {
                    IsValidMatch = (Seq1[P1 - 1] == Seq2[P2 - 1]);
                }
                Data.push_front(typename BaseType::EntryType(Seq1[P1 - 1], Seq2[P2 - 1],
                                                             IsValidMatch));
                P1--;
                P2--;
            }
            // Gap in first sequence
            else if (type == AlignType::_G)
            {
                Data.push_front(
                    typename BaseType::EntryType(Blank, Seq2[P2 - 1], false));
                P2--;
            }
            // Gap in second sequence
            else
            {
                Data.push_front(
                    typename BaseType::EntryType(Seq1[P1 - 1], Blank, false));
                P1--;
            }
        }
    }

    void clearAll()
    {
        if (Matches)
            delete[] Matches;
        Matches = nullptr;

        if (c)
        {
            for (int i = 0; i < (SizeSeq1 + 1) * (SizeSeq2 + 1) + 2; i++)
            {
                delete c[i];
            }
            delete[] c;
        }

        c = nullptr;

        if (hpqueue)
            hpqueue->deleteTable();
        if (hpqueue)
            delete hpqueue;
        hpqueue = nullptr;
    }

public:
    static ScoringSystem getDefaultScoring() { return ScoringSystem(-1, 2, -1); }

    FOGSAASA() : BaseType(getDefaultScoring(), nullptr), Matches(nullptr) {}

    FOGSAASA(ScoringSystem Scoring, MatchFnTy Match = nullptr)
        : BaseType(Scoring, Match), Matches(nullptr) {}

    virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1,
                                                    ContainerType &Seq2)
    {
        AlignedSequence<Ty, Blank> Result;
        // auto t1 = std::chrono::high_resolution_clock::now();
        cacheAllMatches(Seq1, Seq2);
        // auto t2 = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2
        // - t1).count(); std::cout << "caching: " << duration << std::endl;

        // auto t3 = std::chrono::high_resolution_clock::now();
        buildAlignment(Seq1, Seq2, Result);
        // auto t4 = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::microseconds>(t4 -
        // t3).count(); std::cout << "building: " << duration << std::endl;
        clearAll();
        // std::cout << "Expanded nodes: " << expanded << std::endl;
        return Result;
    }

    virtual AlignedSequence<Ty, Blank> getAlignment(std::vector<ContainerType> Seqs)
    {
        assert(Seqs.size() == 2 && "This is a Pairwise Aligner NOT Multiple");
        return getAlignment(Seqs[0], Seqs[1]);
    }
};
