/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2008 Tord Romstad (Glaurung author)
  Copyright (C) 2008-2012 Marco Costalba, Joona Kiiski, Tord Romstad

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(POSITION_H_INCLUDED)
#define POSITION_H_INCLUDED

#include <cassert>

#include "bitboard.h"
#include "types.h"


/// The checkInfo struct is initialized at c'tor time and keeps info used
/// to detect if a move gives check.
class Position;
class Thread;

struct CheckInfo {

  explicit CheckInfo(const Position&);

  Bitboard dcCandidates;
  Bitboard pinned;
  Bitboard checkSq[8];
  Square ksq;
};


/// The StateInfo struct stores information we need to restore a Position
/// object to its previous state when we retract a move. Whenever a move
/// is made on the board (by calling Position::do_move), an StateInfo object
/// must be passed as a parameter.

struct StateInfo {
  Key pawnKey, materialKey;
  Value npMaterial[2];
  int castleRights, rule50, pliesFromNull;
  Score psqScore;
  Square epSquare;

  Key key;
  Bitboard checkersBB;
  PieceType capturedType;
  StateInfo* previous;
};


/// The position data structure. A position consists of the following data:
///
///    * For each piece type, a bitboard representing the squares occupied
///      by pieces of that type.
///    * For each color, a bitboard representing the squares occupied by
///      pieces of that color.
///    * A bitboard of all occupied squares.
///    * A bitboard of all checking pieces.
///    * A 64-entry array of pieces, indexed by the squares of the board.
///    * The current side to move.
///    * Information about the castling rights for both sides.
///    * The initial files of the kings and both pairs of rooks. This is
///      used to implement the Chess960 castling rules.
///    * The en passant square (which is SQ_NONE if no en passant capture is
///      possible).
///    * The squares of the kings for both sides.
///    * Hash keys for the position itself, the current pawn structure, and
///      the current material situation.
///    * Hash keys for all previous positions in the game for detecting
///      repetition draws.
///    * A counter for detecting 50 move rule draws.

class Position {
public:
  Position() {}
  Position(const Position& p) { *this = p; }
  Position(const Position& p, Thread* t) { *this = p; thisThread = t; }
  Position(const std::string& f, bool c960, Thread* t) { from_fen(f, c960, t); }
  void operator=(const Position&);

  // Text input/output
  void from_fen(const std::string& fen, bool isChess960, Thread* th);
  const std::string to_fen() const;
  void print(Move m = MOVE_NONE) const;

  // Position representation
  Bitboard pieces() const;
  Bitboard pieces(Color c) const;
  Bitboard pieces(PieceType pt) const;
  Bitboard pieces(PieceType pt, Color c) const;
  Bitboard pieces(PieceType pt1, PieceType pt2) const;
  Bitboard pieces(PieceType pt1, PieceType pt2, Color c) const;
  Piece piece_on(Square s) const;
  Square king_square(Color c) const;
  Square ep_square() const;
  bool square_empty(Square s) const;
  const Square* piece_list(Color c, PieceType pt) const;
  int piece_count(Color c, PieceType pt) const;

  // Castling
  bool can_castle(CastleRight f) const;
  bool can_castle(Color c) const;
  bool castle_impeded(Color c, CastlingSide s) const;
  Square castle_rook_square(Color c, CastlingSide s) const;

  // Checking
  bool in_check() const;
  Bitboard checkers() const;
  Bitboard discovered_check_candidates() const;
  Bitboard pinned_pieces() const;
  Bitboard all_check_blockers() const;

  // Attacks to/from a given square
  Bitboard attackers_to(Square s) const;
  Bitboard attackers_to(Square s, Bitboard occ) const;
  Bitboard attacks_from(Piece p, Square s) const;
  static Bitboard attacks_from(Piece p, Square s, Bitboard occ);
  template<PieceType> Bitboard attacks_from(Square s) const;
  template<PieceType> Bitboard attacks_from(Square s, Color c) const;

  // Properties of moves
  bool move_gives_check(Move m, const CheckInfo& ci) const;
  bool move_attacks_square(Move m, Square s) const;
  bool pl_move_is_legal(Move m, Bitboard pinned) const;
  bool is_pseudo_legal(const Move m) const;
  bool is_capture(Move m) const;
  bool is_capture_or_promotion(Move m) const;
  bool is_passed_pawn_push(Move m) const;
  Piece piece_moved(Move m) const;
  PieceType captured_piece_type() const;

  // Piece specific
  bool pawn_is_passed(Color c, Square s) const;
  bool pawn_on_7th(Color c) const;
  bool opposite_bishops() const;
  bool bishop_pair(Color c) const;

  // Doing and undoing moves
  void do_move(Move m, StateInfo& st);
  void do_move(Move m, StateInfo& st, const CheckInfo& ci, bool moveIsCheck);
  void undo_move(Move m);
  template<bool Do> void do_null_move(StateInfo& st);

  // Static exchange evaluation
  int see(Move m) const;
  int see_sign(Move m) const;

  // Accessing hash keys
  Key key() const;
  Key exclusion_key() const;
  Key pawn_key() const;
  Key material_key() const;

  // Incremental piece-square evaluation
  Score psq_score() const;
  Score psq_delta(Piece p, Square from, Square to) const;
  Value non_pawn_material(Color c) const;

  // Other properties of the position
  Color side_to_move() const;
  int startpos_ply_counter() const;
  bool is_chess960() const;
  Thread* this_thread() const;
  int64_t nodes_searched() const;
  void set_nodes_searched(int64_t n);
  template<bool SkipRepetition> bool is_draw() const;

  // Position consistency check, for debugging
  bool pos_is_ok(int* failedStep = NULL) const;
  void flip();

  // Global initialization
  static void init();

private:
  // Initialization helpers (used while setting up a position)
  void clear();
  void put_piece(Piece p, Square s);
  void set_castle_right(Color c, Square rfrom);
  bool move_is_legal(const Move m) const;

  // Helper template functions
  template<bool Do> void do_castle_move(Move m);
  template<bool FindPinned, bool EnemiesOnly> Bitboard hidden_checkers() const;

  // Computing hash keys from scratch (for initialization and debugging)
  Key compute_key() const;
  Key compute_pawn_key() const;
  Key compute_material_key() const;

  // Computing incremental evaluation scores and material counts
  Score compute_psq_score() const;
  Value compute_non_pawn_material(Color c) const;

  // Board and pieces
  Piece board[64];             // [square]
  Bitboard byTypeBB[8];        // [pieceType]
  Bitboard byColorBB[2];       // [color]
  int pieceCount[2][8];        // [color][pieceType]
  Square pieceList[2][8][16];  // [color][pieceType][index]
  int index[64];               // [square]

  // Other info
  int castleRightsMask[64];      // [square]
  Square castleRookSquare[2][2]; // [color][side]
  Bitboard castlePath[2][2];     // [color][side]
  StateInfo startState;
  int64_t nodes;
  int startPosPly;
  Color sideToMove;
  Thread* thisThread;
  StateInfo* st;
  int chess960;

  // Static variables
  static Score pieceSquareTable[16][64]; // [piece][square]
  static Key zobrist[2][8][64];          // [color][pieceType][square]/[piece count]
  static Key zobEp[8];                   // [file]
  static Key zobCastle[16];              // [castleRight]
  static Key zobSideToMove;
  static Key zobExclusion;
};

inline int64_t Position::nodes_searched() const {
  return nodes;
}

inline void Position::set_nodes_searched(int64_t n) {
  nodes = n;
}

inline Piece Position::piece_on(Square s) const {
  return board[s];
}

inline Piece Position::piece_moved(Move m) const {
  return board[from_sq(m)];
}

inline bool Position::square_empty(Square s) const {
  return board[s] == NO_PIECE;
}

inline Color Position::side_to_move() const {
  return sideToMove;
}

inline Bitboard Position::pieces() const {
  return byTypeBB[ALL_PIECES];
}

inline Bitboard Position::pieces(Color c) const {
  return byColorBB[c];
}

inline Bitboard Position::pieces(PieceType pt) const {
  return byTypeBB[pt];
}

inline Bitboard Position::pieces(PieceType pt, Color c) const {
  return byTypeBB[pt] & byColorBB[c];
}

inline Bitboard Position::pieces(PieceType pt1, PieceType pt2) const {
  return byTypeBB[pt1] | byTypeBB[pt2];
}

inline Bitboard Position::pieces(PieceType pt1, PieceType pt2, Color c) const {
  return (byTypeBB[pt1] | byTypeBB[pt2]) & byColorBB[c];
}

inline int Position::piece_count(Color c, PieceType pt) const {
  return pieceCount[c][pt];
}

inline const Square* Position::piece_list(Color c, PieceType pt) const {
  return pieceList[c][pt];
}

inline Square Position::ep_square() const {
  return st->epSquare;
}

inline Square Position::king_square(Color c) const {
  return pieceList[c][KING][0];
}

inline bool Position::can_castle(CastleRight f) const {
  return st->castleRights & f;
}

inline bool Position::can_castle(Color c) const {
  return st->castleRights & ((WHITE_OO | WHITE_OOO) << c);
}

inline bool Position::castle_impeded(Color c, CastlingSide s) const {
  return byTypeBB[ALL_PIECES] & castlePath[c][s];
}

inline Square Position::castle_rook_square(Color c, CastlingSide s) const {
  return castleRookSquare[c][s];
}

template<PieceType Pt>
inline Bitboard Position::attacks_from(Square s) const {

  return  Pt == BISHOP || Pt == ROOK ? attacks_bb<Pt>(s, pieces())
        : Pt == QUEEN  ? attacks_from<ROOK>(s) | attacks_from<BISHOP>(s)
        : StepAttacksBB[Pt][s];
}

template<>
inline Bitboard Position::attacks_from<PAWN>(Square s, Color c) const {
  return StepAttacksBB[make_piece(c, PAWN)][s];
}

inline Bitboard Position::attacks_from(Piece p, Square s) const {
  return attacks_from(p, s, byTypeBB[ALL_PIECES]);
}

inline Bitboard Position::attackers_to(Square s) const {
  return attackers_to(s, byTypeBB[ALL_PIECES]);
}

inline Bitboard Position::checkers() const {
  return st->checkersBB;
}

inline bool Position::in_check() const {
  return st->checkersBB != 0;
}

inline Bitboard Position::discovered_check_candidates() const {
  return hidden_checkers<false, true>();
}

inline Bitboard Position::pinned_pieces() const {
  return hidden_checkers<true, true>();
}

inline Bitboard Position::all_check_blockers() const {
  return hidden_checkers<true, false>();
}

inline bool Position::pawn_is_passed(Color c, Square s) const {
  return !(pieces(PAWN, ~c) & passed_pawn_mask(c, s));
}

inline Key Position::key() const {
  return st->key;
}

inline Key Position::exclusion_key() const {
  return st->key ^ zobExclusion;
}

inline Key Position::pawn_key() const {
  return st->pawnKey;
}

inline Key Position::material_key() const {
  return st->materialKey;
}

inline Score Position::psq_delta(Piece p, Square from, Square to) const {
  return pieceSquareTable[p][to] - pieceSquareTable[p][from];
}

inline Score Position::psq_score() const {
  return st->psqScore;
}

inline Value Position::non_pawn_material(Color c) const {
  return st->npMaterial[c];
}

inline bool Position::is_passed_pawn_push(Move m) const {

  return   type_of(piece_moved(m)) == PAWN
        && pawn_is_passed(sideToMove, to_sq(m));
}

inline int Position::startpos_ply_counter() const {
  return startPosPly + st->pliesFromNull; // HACK
}

inline bool Position::opposite_bishops() const {

  return   pieceCount[WHITE][BISHOP] == 1
        && pieceCount[BLACK][BISHOP] == 1
        && opposite_colors(pieceList[WHITE][BISHOP][0], pieceList[BLACK][BISHOP][0]);
}

inline bool Position::bishop_pair(Color c) const {

  return   pieceCount[c][BISHOP] >= 2
        && opposite_colors(pieceList[c][BISHOP][0], pieceList[c][BISHOP][1]);
}

inline bool Position::pawn_on_7th(Color c) const {
  return pieces(PAWN, c) & rank_bb(relative_rank(c, RANK_7));
}

inline bool Position::is_chess960() const {
  return chess960;
}

inline bool Position::is_capture_or_promotion(Move m) const {

  assert(is_ok(m));
  return is_special(m) ? !is_castle(m) : !square_empty(to_sq(m));
}

inline bool Position::is_capture(Move m) const {

  // Note that castle is coded as "king captures the rook"
  assert(is_ok(m));
  return (!square_empty(to_sq(m)) && !is_castle(m)) || is_enpassant(m);
}

inline PieceType Position::captured_piece_type() const {
  return st->capturedType;
}

inline Thread* Position::this_thread() const {
  return thisThread;
}

#endif // !defined(POSITION_H_INCLUDED)
