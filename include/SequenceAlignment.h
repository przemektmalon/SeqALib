//===-- llvm/ADT/SequenceAlignment.h - Sequence Alignment -------*- C++ -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// Provides efficient implementations of different algorithms for sequence
// alignment.
//
//===----------------------------------------------------------------------===//

#ifndef LLVM_ADT_SEQUENCE_ALIGNMENT_H
#define LLVM_ADT_SEQUENCE_ALIGNMENT_H

#include "llvm/ADT/ArrayView.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <limits.h> // INT_MIN
#include <list>
#include <iostream>
#include <fstream>

#define ScoreSystemType int

// Store alignment result here
template <typename Ty, Ty Blank = Ty(0)> class AlignedSequence {
public:
  class Entry {
  private:
    // TODO: change it for a vector<Ty> for Multi-Sequence Alignment
    std::pair<Ty, Ty> Pair;
    bool IsMatchingPair;

  public:
    Entry() { IsMatchingPair = false; }

    Entry(Ty V1, Ty V2) : Pair(V1, V2) { IsMatchingPair = !hasBlank(); }

    Entry(Ty V1, Ty V2, bool Matching)
        : Pair(V1, V2), IsMatchingPair(Matching) {}

    Ty get(size_t index) {
      assert((index == 0 || index == 1) && "Index out of bounds!");
      if (index == 0)
        return Pair.first;
      else
        return Pair.second;
    }

    bool empty() { return (Pair.first == Blank && Pair.second == Blank); }
    bool hasBlank() { return (Pair.first == Blank || Pair.second == Blank); }

    bool match() { return IsMatchingPair; }
    bool mismatch() { return (!IsMatchingPair); }

    Ty getNonBlank() {
      if (Pair.first != Blank)
        return Pair.first;
      else
        return Pair.second;
    }
  };

  std::list<Entry> Data;

  AlignedSequence() {}

  AlignedSequence(const AlignedSequence<Ty, Blank> &Other) : Data(Other.Data) {}
  AlignedSequence(AlignedSequence<Ty, Blank> &&Other)
      : Data(std::move(Other.Data)) {}

  AlignedSequence<Ty> &operator=(const AlignedSequence<Ty, Blank> &Other) {
    Data = Other.Data;
    return (*this);
  }

  void append(const AlignedSequence<Ty, Blank> &Other) {
    Data.insert(Data.end(), Other.Data.begin(), Other.Data.end());
  }

  void splice(AlignedSequence<Ty, Blank> &Other) {
    Data.splice(Data.end(), Other.Data);
  }

  typename std::list<Entry>::iterator begin() { return Data.begin(); }
  typename std::list<Entry>::iterator end() { return Data.end(); }
};

class ScoringSystem {
  ScoreSystemType Gap;
  ScoreSystemType Match;
  ScoreSystemType Mismatch;
  ScoreSystemType GapOpen;
  ScoreSystemType GapExtend;
  bool AllowMismatch;

public:
  ScoringSystem(ScoreSystemType Gap, ScoreSystemType Match) {
    this->Gap = Gap;
    this->Match = Match;
    this->Mismatch = std::numeric_limits<ScoreSystemType>::min();
    this->AllowMismatch = false;
  }

  ScoringSystem(ScoreSystemType Gap, ScoreSystemType Match,
                ScoreSystemType Mismatch, bool AllowMismatch = true) {
    this->Gap = Gap;
    this->Match = Match;
    this->Mismatch = Mismatch;
    this->AllowMismatch = AllowMismatch;
  }

  ScoringSystem(ScoreSystemType GapOpen, ScoreSystemType GapExtend,
                ScoreSystemType Match, ScoreSystemType Mismatch,
                bool AllowMismatch = true) {
    this->GapOpen = GapOpen;
    this->GapExtend = GapExtend;
    this->Match = Match;
    this->Mismatch = Mismatch;
    this->AllowMismatch = AllowMismatch;
  }

  bool getAllowMismatch() { return AllowMismatch; }

  ScoreSystemType getMismatchPenalty() { return Mismatch; }

  ScoreSystemType getGapPenalty() { return Gap; }

  ScoreSystemType getMatchProfit() { return Match; }

  ScoreSystemType getGapOpenPenalty() { return GapOpen; }

  ScoreSystemType getGapExtendPenalty() { return GapExtend; }
};

template <typename ContainerType,
          typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0),
          typename MatchFnTy = std::function<bool(Ty, Ty)>>
class SequenceAligner {
private:
  ScoringSystem Scoring;
  MatchFnTy Match;

public:
  using EntryType = typename AlignedSequence<Ty, Blank>::Entry;

  SequenceAligner(ScoringSystem Scoring, MatchFnTy Match = nullptr)
      : Scoring(Scoring), Match(Match) {}

  ScoringSystem &getScoring() { return Scoring; }

  bool match(Ty Val1, Ty Val2) { return Match(Val1, Val2); }

  MatchFnTy getMatchOperation() { return Match; }

  Ty getBlank() { return Blank; }

  virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq0,
                                                  ContainerType &Seq1) = 0;

  // Force an alignment to be global
  void forceGlobal(ContainerType &Seq1, ContainerType &Seq2,
                   AlignedSequence<Ty, Blank> &Result, int idx1, int idx2,
                   int endIdx1, int endIdx2) {
    auto &Data = Result.Data;

    // Stick on front of sequences up to alignment
    AlignedSequence<Ty, Blank> front;
    for (int i = 0; i < idx1; i++) {
      front.Data.push_back(EntryType(Seq1[i], Blank, false));
    }
    for (int i = 0; i < idx2; i++) {
      front.Data.push_back(EntryType(Blank, Seq2[i], false));
    }

    front.splice(Result);

    // Stick on end of sequences
    AlignedSequence<Ty, Blank> end;
    for (int i = endIdx1; i < Seq1.size(); i++) {
      end.Data.push_back(EntryType(Seq1[i], Blank, false));
    }
    for (int i = endIdx2; i < Seq2.size(); i++) {
      end.Data.push_back(EntryType(Blank, Seq2[i], false));
    }

    front.splice(end);

    Result.Data.clear();

    Result.splice(front);

  }
};

#include "llvm/ADT/SANeedlemanWunsch.h"
#include "llvm/ADT/SADiagonalWindows.h"
#include "llvm/ADT/SAHirschberg.h"
#include "llvm/ADT/SASmithWaterman.h"
#include "llvm/ADT/SAGlobalGotoh.h"
#include "llvm/ADT/SALocalGotoh.h"
#include "llvm/ADT/SAMyersMiller.h"
#include "llvm/ADT/SAFOGSAA.h"
#include "llvm/ADT/CandidateWord.h"
#include "llvm/ADT/MUM.h"
#include "llvm/ADT/SuffixTree.h"
#include "llvm/ADT/SABLAT.h"
#include "llvm/ADT/SAMummer.h"
#include "llvm/ADT/SAGappedBLAT.h"

#endif
