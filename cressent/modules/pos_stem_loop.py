#!/usr/bin/env python3

# this is a modified version of StemLoop-Finder from 
# Pratt, A.A., Torrance, E.L., Kasun, G.W., Stedman, K.M. and de la Higuera, I., 2021. StemLoop-Finder: a tool for the detection of DNA hairpins with conserved motifs. Microbiology Resource Announcements, 10(26), pp.10-1128.
# https://doi.org/10.1128/mra.00424-21
# consider citing this if using this module
# All rights reserved.

import RNA

class PosStemLoop:
    def __init__(self, name, sequence, position, motif, motif_start, reverse_compliment):
        self.name = name
        self.sequence = sequence
        self.string_sequence = str(sequence) + (str(sequence)[0:50])
        self.stem_start = -1
        self.stem_end = -1
        self.stem_length = 0
        self.loop_length = 0
        self.position = position
        self.motif = motif
        self.motif_start = motif_start
        self.folded = ""
        self.score = 0
        self.reverse_compliment = reverse_compliment


    def __repr__(self):
        return repr(self.name)

    def __lt__(self, other):
        return (self.name, self.motif_start) < (other.name, other.motif_start)

    def __eq__(self ,other):
        return (self.name, self.motif_start) == (other.name, other.motif_start)

    def fold(self):
        folded_fragment = RNA.fold_compound(self.sequence)
        (folded_fragment, fragment_mfe) = folded_fragment.mfe()
        self.folded = folded_fragment
        return self.folded

    def find_stem_loop(self, ideal_stem_length, ideal_loop_length):
        pos = -1
        state = 'no_stem'
        stem_count = 0
        loop_count = 0
        stem_match = 0
        #....(((((.......)))))..... <- a loop
        #...(..(...))))).... <- not a loop
        #..(((((.....))))... <- a loop
        while pos < len(self.folded) - 1:
            pos += 1
            c = self.folded[pos]
            if state == 'no_stem':
                stem_count = 0
                if c == '.' or c == ')':
                    continue
                if c == "(":
                    stem_count = 1
                    state = "stem_start"
            elif state == "stem_start":
                if c == '(':
                    stem_count += 1
                    continue
                elif c == ')':
                    state = 'no_stem'
                elif c == '.':
                    loop_count = 1
                    state = 'loop'
            elif state == 'loop':
                if c == '(':
                    stem_count = 1
                    state = 'stem_start'
                if c == '.':
                    loop_count += 1
                    continue
                if c == ')':
                    stem_match = stem_count - 1
                    state = 'stem_end'
            elif state == 'stem_end':
                if c == ')' and stem_match > 1:
                    stem_match -= 1
                elif c == '.':
                    self.stem_len = stem_count - stem_match
                    self.loop_len = loop_count
                    self.stem_end = self.position + pos + 1
                    #self.stem_start = self.position + (self.stem_end - (2 * self.stem_len + self.loop_len) + 2)
                    self.stem_start = self.position + (pos - (2 * self.stem_len + self.loop_len) + 2)
                    if ((self.loop_len > 6) and (self.stem_len > 4)):
                        self.score = abs(self.stem_len - ideal_stem_length) + abs(self.loop_len - ideal_loop_length)
                        stem_loop = self.sequence[(self.stem_start-self.position - 1):(self.stem_end - self.position)]
                        return stem_loop
                    else:
                        self.score = 99

                    state = 'no_stem'
                elif c == '(':
                    self.stem_len = stem_count - stem_match
                    self.loop_len = loop_count
                    self.stem_end = self.position+ pos + 1
                    self.stem_start = self.position + (pos - (2 * self.stem_len + self.loop_len) + 2)
                    if ((self.loop_len > 6) and (self.stem_len > 4)):
                        self.score = abs(self.stem_len - ideal_stem_length) + abs(self.loop_len - ideal_loop_length)
                        stem_loop = self.sequence[(self.stem_start-self.position - 1):(self.stem_end - self.position)]
                        return stem_loop
                    else:
                        self.score = 99

                    state = 'stem_start'
                elif c == ')' and stem_match == 1:
                    self.stem_len = stem_count
                    self.loop_len = loop_count
                    self.stem_end = self.position + pos + 1
                    self.stem_start = self.position + (pos - (2 * self.stem_len + self.loop_len) + 2)
                    if (self.loop_len > 6) and (self.stem_len > 4):
                        self.score = abs(self.stem_len - ideal_stem_length) + abs(self.loop_len - ideal_loop_length)
                        #stem_loop = self.sequence[(self.stem_start-self.position - 1):(self.stem_end - self.position)]
                        stem_loop = self.sequence[(self.stem_start - self.position - 1):(self.stem_end - self.position)]
                        return stem_loop
                    else:
                        self.score = 99
                    state = 'no_stem'
        stem_loop = self.sequence[(self.stem_start-self.position):(self.stem_end - self.position)]
        return stem_loop


    def diff_score_modify(self, ideal_seq):
        if (self.motif == ideal_seq):
            self.score -= -5
            return self.score
        diff_score = 0
        for i in range(len(ideal_seq)):
            if ideal_seq[i] != self.motif[i]:
                diff_score += 1
        self.score += diff_score
        return self.score
