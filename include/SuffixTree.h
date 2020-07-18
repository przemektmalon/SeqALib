#include <iostream>
#include <fstream>
#include <memory>
//<----- Suffix Node ----->
class SuffixNode {
public:
  // Start, end specify the edge by which the node is connected to it's parent
  int start;
  int secondStart;
  std::shared_ptr<int> end;

  SuffixNode **children;  // Children of node
                          // Pointer to pointer as we don't know size yet since
                          // alphabet may be incredibly large
  SuffixNode *suffixLink; // Suffix Link pointer

  int suffixIndex;

  bool leaf;

  // Constructor
  SuffixNode(int newStart, std::shared_ptr<int> newEnd, SuffixNode *root,
             int childrenSize) {
    start = newStart;
    secondStart = newStart;
    end = newEnd;
    suffixLink = root;
    suffixIndex = -1;
    leaf = false;

    this->children = new SuffixNode *[childrenSize+10];
    for (int i = 0; i < childrenSize+10; i++) {
      this->children[i] = nullptr;
    }
  }

  // Destructor
  ~SuffixNode() { delete[] children; }

  SuffixNode *operator=(const SuffixNode &b) {
    this->children = b.children;
    this->start = b.start;
    this->end = b.end;
    this->suffixLink = b.suffixLink;
    this->suffixIndex = b.suffixIndex;
    this->secondStart = b.secondStart;
    this->leaf = b.leaf;
  }
};

//<----- Suffix Tree ----->
template <typename ContainerType,
          typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0),
          typename MatchFnTy = std::function<bool(Ty, Ty)>>
class SuffixTree {
private:
  SuffixNode *lastNewNode = nullptr;
  SuffixNode *activeNode = nullptr;

  int activeEdge = -1;
  int activeLength = 0;

  int remainingSuffixCount = 0;
  std::shared_ptr<int> leafEnd = nullptr;
  // int *rootEnd = nullptr;
  // int *splitEnd = nullptr;
  std::shared_ptr<int> splitEnd = nullptr;
  std::shared_ptr<int> rootEnd = nullptr;
  int seqSize = -1; // Length of input string
  int size1 = 0;    // Length of 1st string

  SuffixNode *root;

  std::vector<Ty> seq;

  std::vector<Ty> alphabet;

  MatchFnTy match;

  using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

  bool seenHash = false;

  bool seenDollar = false;

  bool general = true;

public:
  SuffixTree() {}

  // Builds the generalised suffix tree
  // Generalised suffix tree used for getting the LCS
  void buildTree(ContainerType &Seq1, ContainerType &Seq2, MatchFnTy matchFn) {
    // seq = Seq1 + "#" + Seq2 + "$";
    // seq = Seq1;
    // basic file operations
    match = matchFn;

    for (int i = 0; i < Seq1.size(); i++) {
      seq.push_back(Seq1[i]);
    }

    Ty hash = nullptr;
    Ty dollar = nullptr;
    seq.push_back(hash);

    for (int i = 0; i < Seq2.size(); i++) {
      seq.push_back(Seq2[i]);
    }

    seq.push_back(dollar);

    seqSize = seq.size();
    size1 = Seq1.size() + 1;

    int *a = new int;
    std::shared_ptr<int> newPtr(a);
    rootEnd = newPtr;
    *rootEnd = -1;

    int *b = new int;
    std::shared_ptr<int> newPtrb(b);
    leafEnd = newPtrb;
    *leafEnd = -1;
    root = new SuffixNode(-1, rootEnd, nullptr, seqSize);
    root->suffixLink = root;

    activeNode = root;

    general = true;

    // Begin construction
    /*std::ofstream newfile;
    std::string newPath = "/home/sean/treeChecker.txt";
    newfile.open(newPath, std::ios_base::app);
    newfile << "Attempting.\n";
    newfile << "Seq1 Size: " << Seq1.size() << "\n";
    newfile << "Seq2 Size: " << Seq2.size() << "\n";
    newfile.close();*/
    for (int i = 0; i < seqSize; i++) {
      //newfile.open(newPath, std::ios_base::app);
      //newfile << "Building: " << i << "\n";
      //newfile.close();
      extendSuffixTree(i);
      //newfile.open(newPath, std::ios_base::app);
      //newfile << "Built: " << i << "\n";
      //newfile.close();
    }
    //newfile.open(newPath, std::ios_base::app);
    //newfile << "alphabet size: " << alphabet.size() << "\n";
    //newfile << "is root null? " << (root == nullptr) << "\n";
    //newfile << "We get past extension \n";
    //newfile << "We gonna do set by DFS \n";
    //newfile.close();
    // print(root);
    // newfile.open(newPath, std::ios_base::app);
    // newfile << "We get past printing \n";
    int labelHeight = 0;
    // newfile << "\n";
    // newfile.close();

    setSuffixIndexByDFS(root, labelHeight);
    //newfile.open(newPath, std::ios_base::app);
    //newfile << "Set by DFS \n";
    //newfile.close();
  }

  // Builds the suffix tree
  void buildTree(ContainerType &Seq1, MatchFnTy matchFn) {
    // seq = Seq1;
    for (int i = 0; i < Seq1.size(); i++) {
      seq.push_back(Seq1[i]);
    }
    Ty dollar = nullptr;
    seq.push_back(dollar);

    seqSize = seq.size();
    size1 = seqSize;

    int *a = new int;
    std::shared_ptr<int> newPtr(a);
    rootEnd = newPtr;
    *rootEnd = -1;

    int *b = new int;
    std::shared_ptr<int> newPtrb(b);
    leafEnd = newPtrb;
    *leafEnd = -1;

    root = new SuffixNode(-1, rootEnd, nullptr, seqSize);
    root->suffixLink = root;

    match = matchFn;

    activeNode = root;

    general = false;

    // Begin construction
    for (int i = 0; i < seqSize; i++) {
      extendSuffixTree(i);
    }

    // print(root);
    int labelHeight = 0;
    setSuffixIndexByDFS(root, labelHeight);
  }

  void print(SuffixNode *node) {

    for (int i = 0; i < alphabet.size(); i++) {
      if (node->children[i] != nullptr) {
        // std::string subst = seq.substr(node->children[i]->start,
        //*(node->children[i]->end) -
        // node->children[i]->start + 1);
        // std::cout << seq[node->children[i]->start] << std::endl;
        // std::cout << subst << std::endl;
        // newfile << "start: " << node->children[i]->start << "\n";
        // newfile << "end: " << *(node->children[i]->end) << "\n";
        // newfile << "length: "
        //<< *(node->children[i]->end) - node->children[i]->start + 1
        //<< "\n"
      }
    }
    for (int i = 0; i < alphabet.size(); i++) {
      if (node->children[i] != nullptr) {
        print(node->children[i]);
      }
    }
  }

  void setSuffixIndexByDFS(SuffixNode *n, int labelHeight) {

    /*std::ofstream newfile;
    std::string newPath = "/home/sean/treeChecker.txt";
    newfile.open(newPath, std::ios_base::app);
    newfile << "Setting suffix index\n";
    newfile.close();*/

    if (n == nullptr) {
      return;
    }

    int leaf = 1;
    for (int i = 0; i < alphabet.size(); i++) {
      if (n->children[i] != nullptr) {
        leaf = 0;
        setSuffixIndexByDFS(n->children[i],
                            labelHeight + edgeLength(n->children[i]));
      }
    }

    //newfile.open(newPath, std::ios_base::app);
    //newfile << "Finished Calling Setting\n";
    //newfile.close();

    if (leaf == 1) {

      for (int i = n->start; i <= *(n->end); i++) {
        if (myMatch(seq[i], nullptr, i)) {
          //newfile.open(newPath, std::ios_base::app);
          //newfile << "This nonsense\n";
          //newfile.close();
          n->end = std::shared_ptr<int>(new int);
          *(n->end) = i;
          //newfile.open(newPath, std::ios_base::app);
          //newfile << "Not this nonsense\n";
          //newfile.close();
        }
      }
      n->suffixIndex = seqSize - labelHeight;
      n->leaf = true;
    }
  }

  int edgeLength(SuffixNode *n) {
    if (n == root) {
      return 0;
    }
    return *(n->end) - (n->start) + 1;
  }

  int walkDown(SuffixNode *currNode) {
    // If activeLength is greater than current edge length,
    // set next internal node as activeNode and adjust activeEdge
    // and activeLength accordingly
    if (activeLength >= edgeLength(currNode)) {
      activeEdge += edgeLength(currNode);
      activeLength -= edgeLength(currNode);
      activeNode = currNode;
      return 1;
    }
    return 0;
  }

  void extendSuffixTree(int pos) {
    /*std::ofstream newfile;
    std::string newPath = "/home/sean/treeChecker.txt";
    newfile.open(newPath, std::ios_base::app);
    newfile << "Remaining suffix count: " << remainingSuffixCount << "\n";
    newfile.close();*/
    if (!leafEnd) {
      int *a = new int;
      std::shared_ptr<int> newPtr(a);
      leafEnd = newPtr;
    }
    *leafEnd = pos + 1;
    remainingSuffixCount++;
    lastNewNode = nullptr;

    // Add all suffixes (yet to be added) one by one
    while (remainingSuffixCount > 0) {

      if (activeLength == 0) {
        activeEdge = pos;
      }

      //newfile.open(newPath, std::ios_base::app);
      //newfile << "Remaining suffix count: " << remainingSuffixCount << "\n";
      //newfile << "activeLength: " << activeLength << "\n";
      //newfile.close();

      //newfile.open(newPath, std::ios_base::app);
      //newfile << "Finding\n";
      //newfile.close();
      auto result1 = find(std::begin(alphabet), std::end(alphabet),
                          seq[pos - activeLength], pos);

      if (result1 == std::end(alphabet)) {
        // std::cout << "Can't find " << seq[pos - activeLength] << " in
        // alphabet, adding now" << std::endl;
        alphabet.push_back(seq[pos - activeLength]);
        // std::cout << "Alphabet size: " << alphabet.size() << std::endl;
        activeEdge = alphabet.size() - 1;
        //newfile.open(newPath, std::ios_base::app);
        //newfile << "Not Found: " << activeEdge << "\n";
        //newfile.close();
      } else {
        // std::cout << "Found " << seq[pos - activeLength] << " in alphabet" <<
        // std::endl;
        // Get index of element from iterator
        activeEdge = std::distance(alphabet.begin(), result1);
        //newfile.open(newPath, std::ios_base::app);
        //newfile << "Found: " << activeEdge << "\n";
        //newfile.close();
      }

      // There is no outgoing edge starting with activeEdge from activeNode

      //newfile.open(newPath, std::ios_base::app);
      //newfile << "Is activeNode->children[activeEdge] the culprit?\n";
      //newfile << "activeEdge: " << activeEdge << "\n";
      //newfile.close();
      if (activeNode->children[activeEdge] == nullptr) {
        //newfile.open(newPath, std::ios_base::app);
        //newfile << "Creating new node\n";
        //newfile.close();
        std::shared_ptr<int> sharedLeaf = leafEnd;
        activeNode->children[activeEdge] =
            new SuffixNode(pos, sharedLeaf, root, seqSize);

        //newfile.open(newPath, std::ios_base::app);
        //newfile << "created node \n";
        //newfile.close();

        if (lastNewNode != nullptr) {
          // std::cout << "LAST NEW NODE SUFFIX LINK IS SET" << std::endl;
          lastNewNode->suffixLink = activeNode;
          lastNewNode = nullptr;
        }

        // newfile.open(newPath, std::ios_base::app);
        // newfile << "Created\n";
        // newfile.close();
      }
      // There is an outgoing edge starting with activeEdge from activeNode
      else {
        // std::cout << "Found " << alphabet[activeEdge] << " in activeNode" <<
        // std::endl;
        // newfile.open(newPath, std::ios_base::app);
        // newfile << "Is we good?\n";
        // newfile.close();
        // std::cout << "there is an outgoing edge" << std::endl;
        // Get the next node at the end of edge starting with activeEdge
        SuffixNode *next = activeNode->children[activeEdge];
        // SuffixNode *next = new SuffixNode(activeEdge, &leafEnd, root,
        // seqSize); *next = *activeNode->children[activeEdge]; std::cout <<
        // "***** important ********" << std::endl; std::cout << *next->end << "
        // "
        // << *activeNode->children[activeEdge]->end << std::endl; std::cout <<
        // "activeLength:   << activeLength << "\n";edge is when initialising next: " << activeEdge << std::endl;
        if (walkDown(next)) {
          // Start from next node (the new activeNode)
          // std::cout << "(" << activeNode->start << "," << activeEdge << ","
          // << activeLength << ")" << std::endl;
          continue;
        }
        // newfile.open(newPath, std::ios_base::app);
        // newfile << "Not walkdown?\n";
        // newfile.close();

        // std::cout << "no need for walkdown" << std::endl;
        // Extension rule 3
        // Current character being processed is already on the edge
        // std::cout << "next->start: " << next->start << " activeLength: " <<
        // activeLength << " pos: " << pos << std::endl; std::cout <<
        // "Attempting to match: " << seq[next->start + activeLength] << " and "
        // << seq[pos]
        // << std::endl;
        /*newfile.open(newPath, std::ios_base::app);
        newfile << "attempting to match\n";
        newfile << "next->start: " << next->start << "\n";
        newfile << "activeLength: "  << activeLength << "\n";
        newfile << "pos: " << pos << "\n";
        newfile.close();*/

        // newfile.open(newPath, std::ios_base::app);
        // newfile << "Boutta match?\n";
        // newfile.close();
        if (myMatch(seq[next->start + activeLength], seq[pos], pos)) {
          // newfile.open(newPath, std::ios_base::app);
          // newfile << "Yes match\n";
          // newfile.close();
          // std::cout << "Rule 3" << std::endl;
          //newfile.open(newPath, std::ios_base::app);
          //newfile << "Rule 3\n";
          //newfile.close();
          // If a newly created node waiting for it's
          // suffix link to be set, then set suffix link
          // to current activeLength:   << activeLength << "\n";node
          if (lastNewNode != nullptr && activeNode != root) {
            lastNewNode->suffixLink = activeNode;
            lastNewNode = nullptr;
          }

          activeLength++;

          // Stop all processing in this phase
          // Move onto next phase
          // std::cout << "(" << activeNode->start << "," << activeEdge << ","
          // << activeLength << ")" << std::endl;

          //newfile.open(newPath, std::ios_base::app);
          //newfile << "Break after match\n";
          //newfile.close();

          break;
        }

        //newfile.open(newPath, std::ios_base::app);
        //newfile << "No match\n";
        //newfile.close();

        // Extension rule 2
        // activePoint is in the middle of the edge being traversed
        // and current character being processed is not on the edge
        // So add a new node and a new leaf edge out of that node
        // newfile.open(newPath, std::ios_base::app);
        // newfile << "rule 2\n";
        // newfile.close();
        auto result1 =
            find(std::begin(alphabet), std::end(alphabet), seq[pos], pos);
        int leafEdge = activeEdge;

        if (result1 == std::end(alphabet)) {
          // std::cout << "Can't find " << seq[pos] << " in alphabet, adding
          // now" << std::endl;
          alphabet.push_back(seq[pos]);
          // std::cout << "Alphabet size: " << alphabet.size() << std::endl;
          // activeEdge = alphabet.size();
          leafEdge = alphabet.size() - 1;
          // newfile.open(newPath, std::ios_base::app);
          // newfile << "Not Found\n";
          // newfile.close();
        } else {
          // std::cout << "Found " << seq[pos] << " in alphabet" << std::endl;
          // Get index of element from iterator
          // activeEdge = std::distance(alphabet.begin(), result1);
          leafEdge = std::distance(alphabet.begin(), result1);
          // newfile.open(newPath, std::ios_base::app);
          // newfile << "Found\n";
          // newfile.close();
        }

        // std::cout << "Rule 2" << std::endl;
        std::shared_ptr<int> newTing(new int);
        splitEnd = newTing;
        *splitEnd = next->start + activeLength - 1;

        //newfile.open(newPath, std::ios_base::app);
        //newfile << "newting\n";
        //newfile.close();

        // New internal node
        // std::cout << "creating leaf out of new internal node when activeLength:   << activeLength << "\n";edge
        // is: " << activeEdge << std::endl;

        SuffixNode *split =
            new SuffixNode(next->start, splitEnd, root, seqSize);
        split->secondStart = pos;
        activeNode->children[activeEdge] = split;

        //newfile.open(newPath, std::ios_base::app);
        //newfile << "Split\n";
        //newfile.close();

        // New leaf coming out of new internal node
        std::shared_ptr<int> sharedLeaf = leafEnd;
        split->children[leafEdge] =
            new SuffixNode(pos, sharedLeaf, root, seqSize);

        //newfile.open(newPath, std::ios_base::app);
        //newfile << "leaf\n";
        //newfile.close();

        // New edge coming out of internal node
        next->start += activeLength;
        // std::cout << "NEXT START IS: " << next->start << std::endl;
        auto result2 = find(std::begin(alphabet), std::end(alphabet),
                            seq[next->start], pos);
        int continuingEdge = activeEdge;
        if (result2 == std::end(alphabet)) {
          alphabet.push_back(seq[next->start]);
          continuingEdge = alphabet.size() - 1;
        } else {
          // Get index of element from iterator
          continuingEdge = std::distance(alphabet.begin(), result2);
        }

        //newfile.open(newPath, std::ios_base::app);
        //newfile << "found again\n";
        //newfile.close();
        // std::cout << "CONTINUING EDGE IS: " << continuingEdge << std::endl;

        //newfile.open(newPath, std::ios_base::app);
        //newfile << "The first fuckery\n";
        //newfile.close();

        /*newfile.open(newPath, std::ios_base::app);
        newfile << "Continuing edge: " << continuingEdge << "\n";
        newfile << "Active edge: " << activeEdge << "\n";
        newfile << "Leaf edge: " << leafEdge << "\n";
        newfile.close();*/

        split->children[continuingEdge] =
            new SuffixNode(pos, sharedLeaf, root, seqSize);
        delete split->children[continuingEdge];
        split->children[continuingEdge] = next;

        //newfile.open(newPath, std::ios_base::app);
        //newfile << "The fuckery\n";
        //newfile.close();

        // We have a new internal node
        // If there is any internal node created in last extensions
        // which are still waiting for it's suffix link reset then reset
        if (lastNewNode != nullptr) {
          // suffixLink of lastNewNode points to current newly created node
          lastNewNode->suffixLink = split;
        }

        lastNewNode = split;
        //newfile.open(newPath, std::ios_base::app);
        //newfile << "End of rule 2\n";
        //newfile.close();
        // std::cout << "LAST NEW NODE IS SET AS SPLIT" << std::endl;
      }

      // One suffix got added in tree
      // Decrement the count of suffixes yet to be added
      remainingSuffixCount--;
      if (activeNode == root && activeLength > 0) {
        // std::cout << "activeLength:   << activeLength << "\n";node detected as root" << std::endl;
        activeLength--;
        activeEdge = pos - remainingSuffixCount + 1;

        /*auto result1 = find(std::begin(alphabet), std::end(alphabet),
                            seq[pos - activeLength]);

        if (result1 == std::end(alphabet)) {
          // std::cout << "Can't find " << seq[pos - activeLength] << " in
          // alphabet, adding now" << std::endl;
          alphabet.push_back(seq[pos - activeLength]);
          // std::cout << "Alphabet size: " << alphabet.size() << std::endl;
          activeEdge = alphabet.size() - 1;
        } else {
          // std::cout << "Found " << seq[pos - activeLength] << " in alphabet"
          // << std::endl;
          // Get index of element from iterator
          activeEdge = std::distance(alphabet.begin(), result1);
        }*/
      } else if (activeNode != root) {
        // std::cout << "activeLength:   << activeLength << "\n";node not detected as root" << std::endl;
        activeNode = activeNode->suffixLink;
      }

      // newfile.open(newPath, std::ios_base::app);
      // newfile << "End of it all\n";
      // newfile.close();
      // std::cout << "(" << activeNode->start << "," << activeEdge << "," <<
      // activeLength << ")" << std::endl;
    }
  }

  //<----- Funcs used for LCS ----->
  AlignedSequence<Ty, Blank> getLCS(ContainerType &Seq1, ContainerType &Seq2,
                                    MatchFnTy matchFn) {

    buildTree(Seq1, Seq2, matchFn);
    ContainerType ret;
    int maxHeight = 0;
    int substringStartIndex = 0;
    int secondIndex = 0;

    traverse(root, 0, &maxHeight, &substringStartIndex, &secondIndex);

    AlignedSequence<Ty, Blank> Result;

    for (int i = 0; i < maxHeight; i++) {
      Result.Data.push_back(typename BaseType::EntryType(
          Seq1[substringStartIndex + i], Seq2[secondIndex + i], true));
    }

    return Result;
  }

  void deleteTree() {
    /*std::ofstream newfile;
    std::string newPath = "/home/sean/treeChecker.txt";
    newfile.open(newPath, std::ios_base::app);
    newfile << "Deleting\n";
    newfile.close();*/

    freeSuffixTreeByPostOrder(root);
    // newfile.open(newPath, std::ios_base::app);
    // newfile << "Deleted \n";
    // newfile.close();
    alphabet.clear();
    seq.clear();
  }

  void freeSuffixTreeByPostOrder(SuffixNode *n) {

    if (n == nullptr)
      return;

    for (int i = 0; i < alphabet.size(); i++) {
      if (n->children[i] != nullptr) {
        freeSuffixTreeByPostOrder(n->children[i]);
      }
    }

    delete n;
  }

  int traverse(SuffixNode *n, int labelHeight, int *maxHeight,
               int *substringStartIndex, int *secondIndex) {
    if (n == nullptr) {
      return 0;
    }

    int i = 0;
    int ret = -1;

    if (n->suffixIndex < 0) // If it is an internal node
    {
      for (i = 0; i < alphabet.size(); i++) {
        if (n->children[i] != nullptr) {
          ret =
              traverse(n->children[i], labelHeight + edgeLength(n->children[i]),
                       maxHeight, substringStartIndex, secondIndex);
          if (n->suffixIndex == -1) {
            n->suffixIndex = ret;
          } else if ((n->suffixIndex == -2 && ret == -3) ||
                     (n->suffixIndex == -3 && ret == -2) ||
                     n->suffixIndex == -4) {
            n->suffixIndex = -4; // Mark node as XY
            if (*maxHeight < labelHeight) {
              *maxHeight = labelHeight;
              *substringStartIndex = *(n->end) - labelHeight + 1;
              *secondIndex = n->secondStart - labelHeight - size1;
            }
          }
        }
      }
    } else if (n->suffixIndex > -1 && n->suffixIndex < size1)
      return -2; // Mark node as X
    else if (n->suffixIndex >= size1)
      return -3; // Mark node as Y
    return n->suffixIndex;
  }

  //<----- Funcs used for pattern matching ----->
  template <typename ArrayType>
  int traverseEdge(ArrayType pattern, int idx, int start, int end) {
    int k = 0;
    // Traverse the edge with character by character matching
    for (k = start; k <= end && idx < pattern.size(); k++, idx++) {

      if (!myMatch(seq[k], pattern[idx], 0))
        return -1; // no match
    }

    if (idx == pattern.size()) {

      return 1; // match
    }
    return 0; // more characters yet to match
  }

  int doTraversalToCountLeaf(SuffixNode *n,
                             std::vector<candidateWord> &candidates) {
    if (n == nullptr)
      return 0;
    if (n->suffixIndex > -1) {
      // std::cout << "Found at position: " << n->suffixIndex << std::endl;
      candidateWord cWord;
      cWord.indexSeq1 = n->suffixIndex;
      candidates.push_back(cWord);
      return 1;
    }
    int count = 0;
    int i = 0;
    for (i = 0; i < alphabet.size(); i++) {
      if (n->children[i] != nullptr) {
        count += doTraversalToCountLeaf(n->children[i], candidates);
      }
    }
    return count;
  }

  int countLeaf(SuffixNode *n, std::vector<candidateWord> &candidates) {
    if (n == nullptr)
      return 0;
    return doTraversalToCountLeaf(n, candidates);
  }

  template <typename ArrayType>
  void doTraversal(SuffixNode *n, ArrayType pattern, int idx,
                   std::vector<candidateWord> &candidates) {
    if (n == nullptr) {
      return; // no built tree so return no match
    }
    int res = -1;
    // If node n is not root node, then traverse edge
    // from node n's parent to node n.
    if (n->start != -1) {
      res = traverseEdge(pattern, idx, n->start, *(n->end));
      if (res == -1) // no match
        return;
      if (res == 1) // match
      {
        if (n->suffixIndex > -1) {
          // std::cout << "susbtring count: 1 at position: " << n->suffixIndex
          // << std::endl;
          candidateWord cWord;
          cWord.indexSeq1 = n->suffixIndex;
          candidates.push_back(cWord);
        }

        else
          countLeaf(n, candidates);
        // std::cout << "substring count: " << countLeaf(n, candidates) <<
        // std::endl;
        return;
      }
    }
    // Get the character index to search
    idx = idx + edgeLength(n);
    // If there is an edge from node n going out
    // with current character str[idx], traverse that edge
    auto result1 =
        find(std::begin(alphabet), std::end(alphabet), pattern[idx], *(n->end));
    int index = idx;
    if (result1 != std::end(alphabet)) {
      index = std::distance(alphabet.begin(), result1);
    }
    // std::cout << "index: " << index << std::endl;
    if (n->children[index] != nullptr)
      doTraversal(n->children[index], pattern, idx, candidates);
    else
      return; // no match
  }

  template <typename ArrayType>
  std::vector<candidateWord> getCandidates(ArrayType pattern) {
    std::vector<candidateWord> candidates;
    doTraversal(root, pattern, 0, candidates);
    return candidates;
  }

  bool checkForSubString(ContainerType pattern) {
    int res = doTraversal(root, pattern, 0);
    /*if (res == 1)
        std::cout << "Pattern is a susbtring" << std::endl;
    else
        std::cout << "Pattern is NOT a substring" << std::endl;*/
    return res;
  }

  //<----- Funcs for getting MUMs ----->
  std::vector<MUM> getMUMs(ContainerType &Seq1, ContainerType &Seq2,
                           MatchFnTy matchFn) {

    /*std::ofstream newfile;
    std::string newPath = "/home/sean/treeChecker.txt";
    newfile.open(newPath, std::ios_base::app);
    newfile << "GBUILDING TREES \n";
    newfile.close();*/

    // Build tree
    buildTree(Seq1, Seq2, matchFn);

    //std::ofstream newfile;
    //std::string newPath = "/home/sean/treeChecker.txt";
    //newfile.open(newPath, std::ios_base::app);
    //newfile << "BUILT TREE GETTING MUMS\n";
    //newfile.close();

    // Mark nodes
    std::vector<SuffixNode *> markedNodes;
    std::vector<int> mumStarts;
    std::vector<bool> isSuffix;
    unsigned int length = 0;
    markNodes(root, markedNodes, root, mumStarts, length);
    //newfile.open(newPath, std::ios_base::app);
    //newfile << "Marked Nodes\n";
    //newfile.close();

    // Unmark nodes with suffix links to marked nodes
    for (int i = 0; i < markedNodes.size(); i++) {
      isSuffix.push_back(false);
      for (int j = 0; j < markedNodes.size(); j++) {
        // Skip itself
        if (i != j) {
          // Remove nodes with suffix links to marked nodes
          if (markedNodes[i] == markedNodes[j]->suffixLink) {
            isSuffix[i] = true;
            break;
          }
        }
      }
    }

    //newfile.open(newPath, std::ios_base::app);
    //newfile << "Unmarked Nodes\n";
    //newfile.close();

    std::vector<SuffixNode *> mumNodes;
    std::vector<int> newMumStarts;
    std::vector<MUM> ret;
    for (int i = 0; i < markedNodes.size(); i++) {
      if (!isSuffix[i]) {
        mumNodes.push_back(markedNodes[i]);
        newMumStarts.push_back(mumStarts[i]);
        MUM mum = MUM(mumStarts[i],
                      (markedNodes[i]->secondStart -
                       (*(markedNodes[i]->end) - mumStarts[i] + 1) - size1),
                      (*(markedNodes[i]->end) - mumStarts[i]) + 1);
        ret.push_back(mum);
      }
    }

    //newfile.open(newPath, std::ios_base::app);
    //newfile << "Saved MUMs size: " << ret.size() << "\n";
    //newfile.close();

    // newfile << "SAVED MUMS: size: " << ret.size() << "\n";

    // If no MUM available then simply return the longest common substring
    if (ret.size() == 0) {
      int maxHeight = 0;
      int substringStartIndex = 0;
      int secondIndex = 0;

      traverse(root, 0, &maxHeight, &substringStartIndex, &secondIndex);
      MUM mum = MUM(substringStartIndex, secondIndex, maxHeight);
      ret.push_back(mum);
      /*newfile << "HAD TO GET LCS\n";

      newfile << "idx seq1: " << mum.indexSeq1 << "\n";
      newfile << "idx seq2: " << mum.indexSeq2 << "\n";
      newfile << "mum length: " << mum.length << "\n";
      newfile << "seq1 size: " << Seq1.size() << "\n";
      newfile << "seq2 size: " << Seq2.size() << "\n";*/
    }

    // If any are bad then remove
    for (int i = 0; i < ret.size(); i++) {
      if (ret[i].indexSeq1 < 0 || ret[i].indexSeq2 < 0 ||
          ret[i].indexSeq1 + ret[i].length > Seq1.size() ||
          ret[i].indexSeq2 + ret[i].length > Seq2.size()) {

        /*newfile << "idx seq1: " << ret[i].indexSeq1 << "\n";
        newfile << "idx seq2: " << ret[i].indexSeq2 << "\n";
        newfile << "mum length: " << ret[i].length << "\n";
        newfile << "seq1 size: " << Seq1.size() << "\n";
        newfile << "seq2 size: " << Seq2.size() << "\n";*/
        ret.erase(ret.begin() + i);
        // newfile << "REMOVED\n";
        i--;
      }
    }

    // newfile << "mums size: " << ret.size() << "\n";
    return ret;
  }

  // Recursively loop over all nodes and mark them if they have two leaf
  // children
  void markNodes(SuffixNode *n, std::vector<SuffixNode *> &markedNodes,
                 SuffixNode *startNode, std::vector<int> &mumStarts,
                 unsigned int length) {
    if (n == nullptr) {
      return;
    }

    int noChildren = 0;
    bool isSeq1 = false;
    bool isSeq2 = false;
    for (int i = 0; i < alphabet.size(); i++) {
      if (n->children[i] != nullptr) {
        if (n->children[i]->leaf) {
          noChildren++;

          if (n->children[i]->suffixIndex >= 0 &&
              n->children[i]->suffixIndex < size1) {
            isSeq1 = true;
          } else if (n->children[i]->suffixIndex >= size1) {
            isSeq2 = true;
          }
        } else {
          if (n == root) {
            markNodes(n->children[i], markedNodes, n->children[i], mumStarts,
                      edgeLength(n->children[i]));
          } else {
            markNodes(n->children[i], markedNodes, startNode, mumStarts,
                      length + edgeLength(n->children[i]));
          }
          noChildren = 3;
        }
      }
    }

    // A POTENTIAL MUM
    if (noChildren == 2 && isSeq1 && isSeq2) {
      markedNodes.push_back(n);
      mumStarts.push_back(n->start + edgeLength(n) - length);
    }
  }

  bool myMatch(Ty a, Ty b, int pos) {
    if (a == Blank || a == nullptr || b == Blank || b == nullptr ||
        a == Ty(0) || b == Ty(0)) {
      if (a == nullptr && b == nullptr) {
        if (pos == size1 - 1) {
          return true;
        }
      }
      return false;
    } else {
      return match(a, b);
    }
  }

  template <class InputIt, class T>
  constexpr InputIt find(InputIt first, InputIt last, const T &value, int pos) {
    if (value == nullptr) {
      if (general) {
        if (!seenHash) {
          seenHash = true;
          return last;
        } else {
          if (pos < seq.size() - 1) {
            return first + size1 - 1;
          } else {
            if (!seenDollar) {
              seenDollar = true;
              return last;
            } else {
              return first + seq.size() - 1;
            }
          }
        }
      } else {
        if (!seenDollar) {
          seenDollar = true;
          return last;
        } else {
          return first + seq.size() - 1;
        }
      }
    }

    for (; first != last; ++first) {
      if (myMatch(*first, value, pos)) {
        return first;
      }
    }

    return last;
  }
};
