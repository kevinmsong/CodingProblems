"""
Kevin Song
BMI 214
Project 1: Alignment Algorithms
Start Date: 9/30/17
Submission Date: 10/12/17

Discussed data structure architecture (briefly) with Ankhit Baghel.
Discussed Needleman-Wunsch and Smith-Waterman with Ravina Jain (working through alignment examples on a whiteboard).

Did not share or receive code with/from other classmates.
"""

import sys
import numpy as np
import copy

# Global variables
sequence_a = ""
sequence_b = ""

global_alignment_chosen = None

dx = 0
ex = 0
dy = 0
ey = 0

seq_a_alphabet_length = 0
seq_b_alphabet_length = 0

seq_a_alphabet = ""
seq_b_alphabet = ""

match_matrix = None     # Sequence A's alphabet is along rows, and sequence B's alphabet is along columns.

# Sequence A is along the i axis (rows), and sequence B is along the j axis (columns)
matrix_m = None         # 2D numpy array that holds score values
matrix_ix = None        # 2D numpy array that holds score values
matrix_iy = None        # 2D numpy array that holds score values

# Pointers are represented as 3D arrays, with each cell containing a boolean to indicate pointers 1-7.
ptr_matrix_m = None     # 3D array is 3 items deep, for representing presence/absence of pointers 1-3.
ptr_matrix_ix = None    # 3D array is 2 items deep, for representing presence/absence of pointers 4-5.
ptr_matrix_iy = None    # 3D array is 2 items deep, for representing presence/absence of pointers 6-7.

epsilon = 0.00001


def read_inputs():

    """
    This function reads in the input file, and assigns initial values to some of the program's global variables
    (e.g., sequence_a, sequence_b, gap penalties, lengths and alphabets for each sequence, and the match matrix).

    Inputs: none
    Returns: none
    """

    line_counter = 1

    input_file = open(sys.argv[1], 'r')

    for line in input_file:
        if line_counter == 1:
            global sequence_a
            sequence_a = line[0:(len(line) - 1)]

        elif line_counter == 2:
            global sequence_b
            sequence_b = line[0:(len(line) - 1)]

        elif line_counter == 3:
            global global_alignment_chosen

            if line[0] == "0":
                global_alignment_chosen = True
            elif line[0] == "1":
                global_alignment_chosen = False

        elif line_counter == 4:
            global dx, ex, dy, ey
            split_strings = line.split()

            dx = float(split_strings[0])
            ex = float(split_strings[1])
            dy = float(split_strings[2])
            ey = float(split_strings[3])

        elif line_counter == 5:
            global seq_a_alphabet_length
            seq_a_alphabet_length = int(line)

        elif line_counter == 6:
            global seq_a_alphabet
            seq_a_alphabet = line[0:(len(line) - 1)]

        elif line_counter == 7:
            global seq_b_alphabet_length
            seq_b_alphabet_length = int(line)

        elif line_counter == 8:
            global seq_b_alphabet
            seq_b_alphabet = line[0:(len(line) - 1)]

        elif line_counter >= 9:
            global match_matrix

            if match_matrix is None:
                edge_length = len(seq_a_alphabet)
                match_matrix = [[0 for x in range(edge_length)] for y in range(edge_length)]

            split_strings = line.split()

            if not split_strings:
                break

            row_num = int(split_strings[0])
            col_num = int(split_strings[1])

            match_score = float(split_strings[4])

            match_matrix[row_num - 1][col_num - 1] = match_score

        line_counter += 1

    input_file.close()


# READ INPUTS


read_inputs()


def instantiate_score_matrices():
    """
    This function instantiates the m, ix, and iy matrices as 2D arrays, with zeroes on the leading edges.

    Inputs: none
    Returns: none
    """

    global matrix_m
    global matrix_ix
    global matrix_iy

    matrix_m = np.empty([len(sequence_a) + 1, len(sequence_b) + 1])
    matrix_m[:] = np.nan
    matrix_m[0, :] = 0
    matrix_m[:, 0] = 0

    matrix_ix = np.empty([len(sequence_a) + 1, len(sequence_b) + 1])
    matrix_ix[:] = np.nan
    matrix_ix[0, :] = 0
    matrix_ix[:, 0] = 0

    matrix_iy = np.empty([len(sequence_a) + 1, len(sequence_b) + 1])
    matrix_iy[:] = np.nan
    matrix_iy[0, :] = 0
    matrix_iy[:, 0] = 0


def instantiate_ptr_matrices():
    """
    This function instantiates the pointer matrices for m, ix, and iy as 3D arrays whose cells contain booleans that
    indicate the presence/absence of specific pointers (pointers 1-7, 3 cells deep for m and 2 cells deep for ix/iy).

    Inputs: none
    Returns: none
    """

    global ptr_matrix_m
    global ptr_matrix_ix
    global ptr_matrix_iy

    ptr_matrix_m = np.empty([len(sequence_a) + 1, len(sequence_b) + 1, 3], dtype=bool)
    ptr_matrix_m[:] = False

    ptr_matrix_ix = np.empty([len(sequence_a) + 1, len(sequence_b) + 1, 2], dtype=bool)
    ptr_matrix_ix[:] = False

    ptr_matrix_iy = np.empty([len(sequence_a) + 1, len(sequence_b) + 1, 2], dtype=bool)
    ptr_matrix_iy[:] = False


def get_match_score(matrix_row_num, matrix_col_num):
    """
    This function will calculate the score of a given letter in sequence A (at a certain row number) aligning with
    a given letter in sequence B (at a certain column number). It is typically called when matrix M is being filled.

    Inputs: row and column indices in matrices
    Returns: score from the already-constructed match matrix (constructed during initial data input)
    """

    seq_a_letter = sequence_a[matrix_row_num - 1]
    seq_b_letter = sequence_b[matrix_col_num - 1]

    match_row_index = seq_a_alphabet.index(seq_a_letter)
    match_col_index = seq_b_alphabet.index(seq_b_letter)

    score = match_matrix[match_row_index][match_col_index]

    return score


def fill_matrices():
    """
    This function will iteratively fill in the score matrices (S, Ix, Iy) in parallel while updating the three
    pointer matrices corresponding to the score matrices.

    If global alignment is chosen during file input, this function will fill in the matrices without regard for the
    signs of S, Ix, and Iy's scores. If local alignment was chosen, then any negative scores will be filled in as
    zeroes.

    Inputs: none
    Returns: none
    """

    global matrix_m
    global matrix_ix
    global matrix_iy
    global ptr_matrix_m
    global ptr_matrix_ix
    global ptr_matrix_iy

    rows = len(sequence_a) + 1
    cols = len(sequence_b) + 1

    for row in range(1, rows):
        for col in range(1, cols):
            # Fill M
            score_1 = float(matrix_m[row - 1, col - 1]) + float(get_match_score(row, col))
            score_2 = float(matrix_ix[row - 1, col - 1]) + float(get_match_score(row, col))
            score_3 = float(matrix_iy[row - 1, col - 1]) + float(get_match_score(row, col))

            max_score_m = float(max(score_1, score_2, score_3))

            if global_alignment_chosen:
                matrix_m[row, col] = max_score_m
            else:
                if max_score_m < 0:
                    matrix_m[row, col] = 0
                else:
                    matrix_m[row, col] = max_score_m

            # Update M ptr matrix

            if abs(score_1 - max_score_m) < epsilon:
                ptr_matrix_m[row, col, 0] = True

            if abs(score_2 - max_score_m) < epsilon:
                ptr_matrix_m[row, col, 1] = True

            if abs(score_3 - max_score_m) < epsilon:
                ptr_matrix_m[row, col, 2] = True

            # Fill Ix

            score_4 = float(matrix_m[row - 1, col]) - dy
            score_5 = float(matrix_ix[row - 1, col]) - ey

            max_score_ix = float(max(score_4, score_5))

            if global_alignment_chosen:
                matrix_ix[row, col] = max_score_ix
            else:
                if max_score_ix < 0:
                    matrix_ix[row, col] = 0
                else:
                    matrix_ix[row, col] = max_score_ix

            # Update Ix ptr matrix

            if abs(score_4 - max_score_ix) < epsilon:
                ptr_matrix_ix[row, col, 0] = True

            if abs(score_5 - max_score_ix) < epsilon:
                ptr_matrix_ix[row, col, 1] = True

            # Fill Iy

            score_6 = float(matrix_m[row, col - 1]) - dx
            score_7 = float(matrix_iy[row, col - 1]) - ex

            max_score_iy = float(max(score_6, score_7))

            if global_alignment_chosen:
                matrix_iy[row, col] = max_score_iy
            else:
                if max_score_iy < 0:
                    matrix_iy[row, col] = 0
                else:
                    matrix_iy[row, col] = max_score_iy

            # Update Iy ptr matrix

            if abs(score_6 - max_score_iy) < epsilon:
                ptr_matrix_iy[row, col, 0] = True

            if abs(score_7 - max_score_iy) < epsilon:
                ptr_matrix_iy[row, col, 1] = True


def get_highest_score_global():
    """
    This function will return the indices and score of the edge-aligned cell in M with the highest score.
    This function is called when global pair-ends-free gap alignment is performed.

    Input: none
    Returns: row_num, col_num, matrix_code: lower-right-hand index and code indicating 0 for M, 1 for Ix, and 2 for Iy.
    """

    final_row_index = len(sequence_a)
    final_col_index = len(sequence_b)

    max_starting_score = -np.inf

    for row in range(0, final_row_index + 1):
        if matrix_m[row, final_col_index] > max_starting_score:
            max_starting_score = matrix_m[row, final_col_index]

    for col in range(0, final_col_index + 1):
        if matrix_m[final_row_index, col] > max_starting_score:
            max_starting_score = matrix_m[final_row_index, col]

    for row in range(0, final_row_index + 1):
        if abs(matrix_m[row, final_col_index] - max_starting_score) < epsilon:
            return row, final_col_index, 0

    for col in range(0, final_col_index + 1):
        if abs(matrix_m[final_row_index, col] - max_starting_score) < epsilon:
            return final_row_index, col, 0


def get_highest_score_local():
    """
    This function will return the indices and score of the cell with the highest score in M. This function is called
    when global alignment is performed.

    Input: none
    Returns: row_num, col_num, matrix_code: highest-scoring cell's indices in M (so matrix_code will always be returned
    as 0).
    """

    highest_score = -np.inf
    row_num = None
    col_num = None

    for row in range(0, len(sequence_a) + 1):
        for col in range(0, len(sequence_b) + 1):
            if matrix_m[row, col] > highest_score:
                highest_score = matrix_m[row, col]
                row_num = row
                col_num = col

    return row_num, col_num, 0


def get_matrix_score(row, col, matrix_code):
    """
    This function will return the score in a given matrix M, Ix, or Iy.

    Inputs: row and column indices of a specific matrix. matrix_code = 0 for M, 1 for Ix, and 2 for Iy.
    Returns: score at the indicated row/column indices of a specific matrix.
    """

    if matrix_code is 0:
        return matrix_m[row, col]
    elif matrix_code is 1:
        return matrix_ix[row, col]
    elif matrix_code is 2:
        return matrix_iy[row, col]
    else:
        return None


def check_pointers(row, col, matrix_code):
    """
    This checks which pointers are present in the pointer matrix corresponding to a specific matrix M, Ix, or Iy.
    Pointers 1, 2, and 3 correspond to the three M matrix pointers. Pointers 4 and 5 correspond to the two Ix pointers,
    and pointers 6 and 7 correspond to the two Iy pointers.

    Inputs: row and column indices of a specific matrix M, Ix, or Iy. matrix_code = 0 for M, 1 for Ix, and 2 for Iy.
    Returns: seven booleans (presence/absence of pointers 1-7)
    """

    pointer_1 = False
    pointer_2 = False
    pointer_3 = False
    pointer_4 = False
    pointer_5 = False
    pointer_6 = False
    pointer_7 = False

    if matrix_code is 0:
        if ptr_matrix_m[row, col, 0]:
            pointer_1 = True
        if ptr_matrix_m[row, col, 1]:
            pointer_2 = True
        if ptr_matrix_m[row, col, 2]:
            pointer_3 = True
    elif matrix_code is 1:
        if ptr_matrix_ix[row, col, 0]:
            pointer_4 = True
        if ptr_matrix_ix[row, col, 1]:
            pointer_5 = True
    elif matrix_code is 2:
        if ptr_matrix_iy[row, col, 0]:
            pointer_6 = True
        if ptr_matrix_iy[row, col, 1]:
            pointer_7 = True

    return pointer_1, pointer_2, pointer_3, pointer_4, pointer_5, pointer_6, pointer_7


def get_sequence_match_score(string_a, string_b):
    """
    This function will calculate the match score for the final output of an alignment of sequence A with sequence B.

    Inputs: substrings of sequence A and sequence B (string_a, string_b)
    Returns: score count for that pair alignment
    """
    score_count = 0

    for str_index in range(0, len(string_a)):
        if string_a[str_index].isalpha() and string_b[str_index].isalpha():
            match_row_index = seq_a_alphabet.index(string_a[str_index])
            match_col_index = seq_b_alphabet.index(string_b[str_index])
            score = match_matrix[match_row_index][match_col_index]
            score_count += score

    for str_index in range(1, len(string_a) - 1):
        if string_a[str_index] is "_" and string_a[str_index - 1] is not "_":
            score_count -= dx
        elif string_a[str_index] is "_" and string_a[str_index - 1] is "_":
            score_count -= ex
        if string_b[str_index] is "_" and string_b[str_index - 1] is not "_":
            score_count -= dy
        elif string_b[str_index] is "_" and string_b[str_index - 1] is "_":
            score_count -= ey

    return score_count


def traceback_rec(row_index, col_index, matrix_code, list_index, list_of_word_pairs):
    """
    This function will perform traceback, and emit letters and/or gaps, one pair at a time.
    It will use recursion, and will store emitted sequence strings in a list of lists, [[string1, string2],
    [string1, string2]]. It will fork when necessary, and keep track of the current row/column indices,
    the matrix that it is currently in (matrix_code 0 = M, 1 = Ix, 2 = Iy). It will also keep track of the index
    of the list that it is currently emitting.

    Inputs: row and column indices of a specific matrix. matrix_code = 0 for M, 1 for Ix, and 2 for Iy. list_index =
    the index of the list within the list of lists that the algorithm is currently emitting. list_of_word_pairs =
    [[string1, string2], [string1, string2], ...] that is currently being built up.
    Returns: a completed list of lists (seq-A:seq-B word-pairs)
    """

    if row_index is 0 and col_index is 0:
        return list_index

    if abs(get_matrix_score(row_index, col_index, matrix_code) - 0) < epsilon:
        return list_index

    else:
        ptr_1, ptr_2, ptr_3, ptr_4, ptr_5, ptr_6, ptr_7 = check_pointers(row_index, col_index, matrix_code)

        word_pair = list(list_of_word_pairs[list_index])

        string_a = copy.copy(word_pair[0])
        string_b = copy.copy(word_pair[1])

        if matrix_code is 0:  # Emit A[i]:B[j]
            if ptr_1 and not ptr_2 and not ptr_3:
                string_a = sequence_a[row_index - 1] + string_a
                string_b = sequence_b[col_index - 1] + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index - 1, col_index - 1, 0, list_index, list_of_word_pairs)

            elif not ptr_1 and ptr_2 and not ptr_3:
                string_a = sequence_a[row_index - 1] + string_a
                string_b = sequence_b[col_index - 1] + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index - 1, col_index - 1, 1, list_index, list_of_word_pairs)

            elif not ptr_1 and not ptr_2 and ptr_3:
                string_a = sequence_a[row_index - 1] + string_a
                string_b = sequence_b[col_index - 1] + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index - 1, col_index - 1, 2, list_index, list_of_word_pairs)

            elif ptr_1 and ptr_2 and not ptr_3:
                ptr_1_word_pair_copy = list(word_pair)
                ptr_2_word_pair_copy = list(word_pair)

                string_a_1 = ptr_1_word_pair_copy[0]
                string_b_1 = ptr_1_word_pair_copy[1]

                string_a_1 = sequence_a[row_index - 1] + string_a_1
                string_b_1 = sequence_b[col_index - 1] + string_b_1

                list_of_word_pairs[list_index][0] = string_a_1
                list_of_word_pairs[list_index][1] = string_b_1

                list_index = traceback_rec(row_index - 1, col_index - 1, 0, list_index, list_of_word_pairs)

                string_a_2 = ptr_2_word_pair_copy[0]
                string_b_2 = ptr_2_word_pair_copy[1]

                string_a_2 = sequence_a[row_index - 1] + string_a_2
                string_b_2 = sequence_b[col_index - 1] + string_b_2

                list_of_word_pairs.insert(list_index + 1, [string_a_2, string_b_2])

                list_index = traceback_rec(row_index - 1, col_index - 1, 1, list_index + 1, list_of_word_pairs)

                return list_index

            elif ptr_1 and not ptr_2 and ptr_3:
                ptr_1_word_pair_copy = list(word_pair)
                ptr_3_word_pair_copy = list(word_pair)

                string_a_1 = ptr_1_word_pair_copy[0]
                string_b_1 = ptr_1_word_pair_copy[1]

                string_a_1 = sequence_a[row_index - 1] + string_a_1
                string_b_1 = sequence_b[col_index - 1] + string_b_1

                list_of_word_pairs[list_index][0] = string_a_1
                list_of_word_pairs[list_index][1] = string_b_1

                list_index = traceback_rec(row_index - 1, col_index - 1, 0, list_index, list_of_word_pairs)

                string_a_3 = ptr_3_word_pair_copy[0]
                string_b_3 = ptr_3_word_pair_copy[1]

                string_a_3 = sequence_a[row_index - 1] + string_a_3
                string_b_3 = sequence_b[col_index - 1] + string_b_3

                list_of_word_pairs.insert(list_index + 1, [string_a_3, string_b_3])

                list_index = traceback_rec(row_index - 1, col_index - 1, 2, list_index + 1, list_of_word_pairs)

                return list_index

            elif not ptr_1 and ptr_2 and ptr_3:
                ptr_2_word_pair_copy = list(word_pair)
                ptr_3_word_pair_copy = list(word_pair)

                string_a_2 = ptr_2_word_pair_copy[0]
                string_b_2 = ptr_2_word_pair_copy[1]

                string_a_2 = sequence_a[row_index - 1] + string_a_2
                string_b_2 = sequence_b[col_index - 1] + string_b_2

                list_of_word_pairs[list_index][0] = string_a_2
                list_of_word_pairs[list_index][1] = string_b_2

                list_index = traceback_rec(row_index - 1, col_index - 1, 1, list_index, list_of_word_pairs)

                string_a_3 = ptr_3_word_pair_copy[0]
                string_b_3 = ptr_3_word_pair_copy[1]

                string_a_3 = sequence_a[row_index - 1] + string_a_3
                string_b_3 = sequence_b[col_index - 1] + string_b_3

                list_of_word_pairs.insert(list_index + 1, [string_a_3, string_b_3])

                list_index = traceback_rec(row_index - 1, col_index - 1, 2, list_index + 1, list_of_word_pairs)

                return list_index

            elif ptr_1 and ptr_2 and ptr_3:
                ptr_1_word_pair_copy = list(word_pair)
                ptr_2_word_pair_copy = list(word_pair)
                ptr_3_word_pair_copy = list(word_pair)

                string_a_1 = ptr_1_word_pair_copy[0]
                string_b_1 = ptr_1_word_pair_copy[1]

                string_a_1 = sequence_a[row_index - 1] + string_a_1
                string_b_1 = sequence_b[col_index - 1] + string_b_1

                list_of_word_pairs[list_index][0] = string_a_1
                list_of_word_pairs[list_index][1] = string_b_1

                list_index = traceback_rec(row_index - 1, col_index - 1, 0, list_index, list_of_word_pairs)

                string_a_2 = ptr_2_word_pair_copy[0]
                string_b_2 = ptr_2_word_pair_copy[1]

                string_a_2 = sequence_a[row_index - 1] + string_a_2
                string_b_2 = sequence_b[col_index - 1] + string_b_2

                list_of_word_pairs.insert(list_index + 1, [string_a_2, string_b_2])

                list_index = traceback_rec(row_index - 1, col_index - 1, 1, list_index + 1, list_of_word_pairs)

                string_a_3 = ptr_3_word_pair_copy[0]
                string_b_3 = ptr_3_word_pair_copy[1]

                string_a_3 = sequence_a[row_index - 1] + string_a_3
                string_b_3 = sequence_b[col_index - 1] + string_b_3

                list_of_word_pairs.insert(list_index + 1, [string_a_3, string_b_3])

                list_index = traceback_rec(row_index - 1, col_index - 1, 2, list_index + 1, list_of_word_pairs)

                return list_index

        elif matrix_code is 1:  # Emit A[i]: "_"
            if ptr_4 and not ptr_5:
                string_a = sequence_a[row_index - 1] + string_a
                string_b = "_" + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index - 1, col_index, 0, list_index, list_of_word_pairs)

            elif not ptr_4 and ptr_5:
                string_a = sequence_a[row_index - 1] + string_a
                string_b = "_" + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index - 1, col_index, 1, list_index, list_of_word_pairs)

            elif ptr_4 and ptr_5:
                ptr_4_word_pair_copy = list(word_pair)
                ptr_5_word_pair_copy = list(word_pair)

                string_a_4 = ptr_4_word_pair_copy[0]
                string_b_4 = ptr_4_word_pair_copy[1]

                string_a_4 = sequence_a[row_index - 1] + string_a_4
                string_b_4 = "_" + string_b_4

                list_of_word_pairs[list_index][0] = string_a_4
                list_of_word_pairs[list_index][1] = string_b_4

                list_index = traceback_rec(row_index - 1, col_index, 0, list_index, list_of_word_pairs)

                string_a_5 = ptr_5_word_pair_copy[0]
                string_b_5 = ptr_5_word_pair_copy[1]

                string_a_5 = sequence_a[row_index - 1] + string_a_5
                string_b_5 = "_" + string_b_5

                list_of_word_pairs.insert(list_index + 1, [string_a_5, string_b_5])

                list_index = traceback_rec(row_index - 1, col_index, 1, list_index + 1, list_of_word_pairs)

                return list_index

        elif matrix_code is 2:  # Emit "_":B[j]
            if ptr_6 and not ptr_7:
                string_a = "_" + string_a
                string_b = sequence_b[col_index - 1] + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index, col_index - 1, 0, list_index, list_of_word_pairs)

            elif not ptr_6 and ptr_7:
                string_a = "_" + string_a
                string_b = sequence_b[col_index - 1] + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index, col_index - 1, 2, list_index, list_of_word_pairs)

            elif ptr_6 and ptr_7:
                ptr_6_word_pair_copy = list(word_pair)
                ptr_7_word_pair_copy = list(word_pair)

                string_a_6 = ptr_6_word_pair_copy[0]
                string_b_6 = ptr_6_word_pair_copy[1]

                string_a_6 = "_" + string_a_6
                string_b_6 = sequence_b[col_index - 1] + string_b_6

                list_of_word_pairs[list_index][0] = string_a_6
                list_of_word_pairs[list_index][1] = string_b_6

                list_index = traceback_rec(row_index, col_index - 1, 0, list_index, list_of_word_pairs)

                string_a_7 = ptr_7_word_pair_copy[0]
                string_b_7 = ptr_7_word_pair_copy[1]

                string_a_7 = "_" + string_a_7
                string_b_7 = sequence_b[col_index - 1] + string_b_7

                list_of_word_pairs.insert(list_index + 1, [string_a_7, string_b_7])

                list_index = traceback_rec(row_index, col_index - 1, 2, list_index + 1, list_of_word_pairs)

                return list_index


def trim_end_gaps(word_pairs):
    """
    This function will remove end gaps for global alignment outputs.

    Inputs: [[list of lists]] of sequence A : sequence B alignment strings
    Output: Modified [[list of lists]], with removed end gaps
    """
    for l in word_pairs:
        while (l[0][0] is "_" and l[1][0] is not "_") or (l[0][0] is not "_" and l[1][0] is "_"):
            l[0] = l[0][1:]
            l[1] = l[1][1:]

        while (l[0][-1] is "_" and l[1][-1] is not "_") or (l[0][-1] is not "_" and l[1][-1] is "_"):
            l[0] = l[0][0:(len(l[0]) - 1)]
            l[1] = l[1][0:(len(l[1]) - 1)]


def write_to_output_file(score, word_list_of_lists):
    """
    This function will write the calculated alignment score to an output file. It will also iterate over a list of
    lists that contains word pairs, and write each word pair to the output file.

    Inputs: score: calculated alignment score for word pairs in the final list of lists,
    word_list_of_lists: final list of lists that contains word pairs aligned by the algorithm
    Returns: none
    """

    output_file = open(sys.argv[2], "w+")

    output_file.write(str(score))

    output_file.write("\n")

    for element_list in word_list_of_lists:
        output_file.write("\n")
        output_file.write(element_list[0])
        output_file.write("\n")
        output_file.write(element_list[1])
        output_file.write("\n")


def run_alignment():
    """
    This function serves as a pseudo-"run()" function for the alignment program. It is responsible for calling other
    functions to perform the alignment.

    Inputs: none
    Returns: none
    """

    word_list_of_lists = [["", ""]]
    list_index = 0

    instantiate_score_matrices()
    instantiate_ptr_matrices()
    fill_matrices()

    if global_alignment_chosen:
        row_num, col_num, matrix_code = get_highest_score_global()
        traceback_rec(row_num, col_num, matrix_code, list_index, word_list_of_lists)
        trim_end_gaps(word_list_of_lists)

    else:
        row_num, col_num, matrix_code = get_highest_score_local()
        traceback_rec(row_num, col_num, matrix_code, list_index, word_list_of_lists)

    word_list_of_lists.sort()

    output = word_list_of_lists
    duplicate_removed_output = [output[i] for i in range(len(output)) if i == 0 or output[i] != output[i - 1]]

    alignment_score = get_sequence_match_score(duplicate_removed_output[0][0], duplicate_removed_output[0][1])

    write_to_output_file(alignment_score, duplicate_removed_output)


# INITIATE PROGRAM

if global_alignment_chosen is True or global_alignment_chosen is False:
    run_alignment()
else:
    print "Error: no alignment method defined"
rd_pair_copy[0]
                string_b_3 = ptr_3_word_pair_copy[1]

                string_a_3 = sequence_a[row_index - 1] + string_a_3
                string_b_3 = sequence_b[col_index - 1] + string_b_3

                list_of_word_pairs.insert(list_index + 1, [string_a_3, string_b_3])

                list_index = traceback_rec(row_index - 1, col_index - 1, 2, list_index + 1, list_of_word_pairs)

                return list_index

            elif not ptr_1 and ptr_2 and ptr_3:
                ptr_2_word_pair_copy = list(word_pair)
                ptr_3_word_pair_copy = list(word_pair)

                string_a_2 = ptr_2_word_pair_copy[0]
                string_b_2 = ptr_2_word_pair_copy[1]

                string_a_2 = sequence_a[row_index - 1] + string_a_2
                string_b_2 = sequence_b[col_index - 1] + string_b_2

                list_of_word_pairs[list_index][0] = string_a_2
                list_of_word_pairs[list_index][1] = string_b_2

                list_index = traceback_rec(row_index - 1, col_index - 1, 1, list_index, list_of_word_pairs)

                string_a_3 = ptr_3_word_pair_copy[0]
                string_b_3 = ptr_3_word_pair_copy[1]

                string_a_3 = sequence_a[row_index - 1] + string_a_3
                string_b_3 = sequence_b[col_index - 1] + string_b_3

                list_of_word_pairs.insert(list_index + 1, [string_a_3, string_b_3])

                list_index = traceback_rec(row_index - 1, col_index - 1, 2, list_index + 1, list_of_word_pairs)

                return list_index

            elif ptr_1 and ptr_2 and ptr_3:
                ptr_1_word_pair_copy = list(word_pair)
                ptr_2_word_pair_copy = list(word_pair)
                ptr_3_word_pair_copy = list(word_pair)

                string_a_1 = ptr_1_word_pair_copy[0]
                string_b_1 = ptr_1_word_pair_copy[1]

                string_a_1 = sequence_a[row_index - 1] + string_a_1
                string_b_1 = sequence_b[col_index - 1] + string_b_1

                list_of_word_pairs[list_index][0] = string_a_1
                list_of_word_pairs[list_index][1] = string_b_1

                list_index = traceback_rec(row_index - 1, col_index - 1, 0, list_index, list_of_word_pairs)

                string_a_2 = ptr_2_word_pair_copy[0]
                string_b_2 = ptr_2_word_pair_copy[1]

                string_a_2 = sequence_a[row_index - 1] + string_a_2
                string_b_2 = sequence_b[col_index - 1] + string_b_2

                list_of_word_pairs.insert(list_index + 1, [string_a_2, string_b_2])

                list_index = traceback_rec(row_index - 1, col_index - 1, 1, list_index + 1, list_of_word_pairs)

                string_a_3 = ptr_3_word_pair_copy[0]
                string_b_3 = ptr_3_word_pair_copy[1]

                string_a_3 = sequence_a[row_index - 1] + string_a_3
                string_b_3 = sequence_b[col_index - 1] + string_b_3

                list_of_word_pairs.insert(list_index + 1, [string_a_3, string_b_3])

                list_index = traceback_rec(row_index - 1, col_index - 1, 2, list_index + 1, list_of_word_pairs)

                return list_index

        elif matrix_code is 1:  # Emit A[i]: "_"
            if ptr_4 and not ptr_5:
                string_a = sequence_a[row_index - 1] + string_a
                string_b = "_" + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index - 1, col_index, 0, list_index, list_of_word_pairs)

            elif not ptr_4 and ptr_5:
                string_a = sequence_a[row_index - 1] + string_a
                string_b = "_" + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index - 1, col_index, 1, list_index, list_of_word_pairs)

            elif ptr_4 and ptr_5:
                ptr_4_word_pair_copy = list(word_pair)
                ptr_5_word_pair_copy = list(word_pair)

                string_a_4 = ptr_4_word_pair_copy[0]
                string_b_4 = ptr_4_word_pair_copy[1]

                string_a_4 = sequence_a[row_index - 1] + string_a_4
                string_b_4 = "_" + string_b_4

                list_of_word_pairs[list_index][0] = string_a_4
                list_of_word_pairs[list_index][1] = string_b_4

                list_index = traceback_rec(row_index - 1, col_index, 0, list_index, list_of_word_pairs)

                string_a_5 = ptr_5_word_pair_copy[0]
                string_b_5 = ptr_5_word_pair_copy[1]

                string_a_5 = sequence_a[row_index - 1] + string_a_5
                string_b_5 = "_" + string_b_5

                list_of_word_pairs.insert(list_index + 1, [string_a_5, string_b_5])

                list_index = traceback_rec(row_index - 1, col_index, 1, list_index + 1, list_of_word_pairs)

                return list_index

        elif matrix_code is 2:  # Emit "_":B[j]
            if ptr_6 and not ptr_7:
                string_a = "_" + string_a
                string_b = sequence_b[col_index - 1] + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index, col_index - 1, 0, list_index, list_of_word_pairs)

            elif not ptr_6 and ptr_7:
                string_a = "_" + string_a
                string_b = sequence_b[col_index - 1] + string_b

                list_of_word_pairs[list_index][0] = string_a
                list_of_word_pairs[list_index][1] = string_b

                return traceback_rec(row_index, col_index - 1, 2, list_index, list_of_word_pairs)

            elif ptr_6 and ptr_7:
                ptr_6_word_pair_copy = list(word_pair)
                ptr_7_word_pair_copy = list(word_pair)

                string_a_6 = ptr_6_word_pair_copy[0]
                string_b_6 = ptr_6_word_pair_copy[1]

                string_a_6 = "_" + string_a_6
                string_b_6 = sequence_b[col_index - 1] + string_b_6

                list_of_word_pairs[list_index][0] = string_a_6
                list_of_word_pairs[list_index][1] = string_b_6

                list_index = traceback_rec(row_index, col_index - 1, 0, list_index, list_of_word_pairs)

                string_a_7 = ptr_7_word_pair_copy[0]
                string_b_7 = ptr_7_word_pair_copy[1]

                string_a_7 = "_" + string_a_7
                string_b_7 = sequence_b[col_index - 1] + string_b_7

                list_of_word_pairs.insert(list_index + 1, [string_a_7, string_b_7])

                list_index = traceback_rec(row_index, col_index - 1, 2, list_index + 1, list_of_word_pairs)

                return list_index


def trim_end_gaps(word_pairs):
    """
    This function will remove end gaps for global alignment outputs.

    Inputs: [[list of lists]] of sequence A : sequence B alignment strings
    Output: Modified [[list of lists]], with removed end gaps
    """
    for l in word_pairs:
        while (l[0][0] is "_" and l[1][0] is not "_") or (l[0][0] is not "_" and l[1][0] is "_"):
            l[0] = l[0][1:]
            l[1] = l[1][1:]

        while (l[0][-1] is "_" and l[1][-1] is not "_") or (l[0][-1] is not "_" and l[1][-1] is "_"):
            l[0] = l[0][0:(len(l[0]) - 1)]
            l[1] = l[1][0:(len(l[1]) - 1)]


def write_to_output_file(score, word_list_of_lists):
    """
    This function will write the calculated alignment score to an output file. It will also iterate over a list of
    lists that contains word pairs, and write each word pair to the output file.

    Inputs: score: calculated alignment score for word pairs in the final list of lists,
    word_list_of_lists: final list of lists that contains word pairs aligned by the algorithm
    Returns: none
    """

    output_file = open(sys.argv[2], "w+")

    output_file.write(str(score))

    output_file.write("\n")

    for element_list in word_list_of_lists:
        output_file.write("\n")
        output_file.write(element_list[0])
        output_file.write("\n")
        output_file.write(element_list[1])
        output_file.write("\n")


def run_alignment():
    """
    This function serves as a pseudo-"run()" function for the alignment program. It is responsible for calling other
    functions to perform the alignment.

    Inputs: none
    Returns: none
    """

    word_list_of_lists = [["", ""]]
    list_index = 0

    instantiate_score_matrices()
    instantiate_ptr_matrices()
    fill_matrices()

    if global_alignment_chosen:
        row_num, col_num, matrix_code = get_highest_score_global()
        traceback_rec(row_num, col_num, matrix_code, list_index, word_list_of_lists)
        trim_end_gaps(word_list_of_lists)

    else:
        row_num, col_num, matrix_code = get_highest_score_local()
        traceback_rec(row_num, col_num, matrix_code, list_index, word_list_of_lists)

    word_list_of_lists.sort()

    output = word_list_of_lists
    duplicate_removed_output = [output[i] for i in range(len(output)) if i == 0 or output[i] != output[i - 1]]

    alignment_score = get_sequence_match_score(duplicate_removed_output[0][0], duplicate_removed_output[0][1])

    write_to_output_file(alignment_score, duplicate_removed_output)


# INITIATE PROGRAM

if global_alignment_chosen is True or global_alignment_chosen is False:
    run_alignment()
else:
    print "Error: no alignment method defined"
