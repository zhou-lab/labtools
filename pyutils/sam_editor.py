import pysam
from sequence_operations_common import reverse_complement

class sam_editor:

    def __init__(self, sam_obj):
        self.sam_obj = sam_obj
        if self.sam_obj.flag & 16 and self.sam_obj.query_qualities is not None:
            q = self.sam_obj.query_qualities[:]
            self.sam_obj.query_sequence = reverse_complement(self.sam_obj.query_sequence)
            self.sam_obj.query_qualities = q[::-1]

    def getHead(self):
        return self.sam_obj.query_name

    def getSeq(self):
        return self.sam_obj.query_sequence

    def getQual(self):
        return self.sam_obj.query_qualities

    def trimming(self, start, end):
        if start > end or end > self.sam_obj.query_length:
            raise Exception("The numbers you entered don't make sense. Come up with more informative error message.")
        self.trimming_methyl_tuple(start, end)
        q = self.sam_obj.query_qualities
        self.sam_obj.query_sequence = self.sam_obj.query_sequence[start:end]
        self.sam_obj.query_qualities = q[start:end]
        # self.sam_obj.query_length = end - start

    def unstring_methyl_tuple(self):
        tuples = self.sam_obj.get_tags()
        MM_tag = ''
        ML_tag = ''
        for tuple_single in tuples:
            if tuple_single[0] == 'MM':
                MM_tag = tuple_single[1]
            elif tuple_single[0] == 'ML':
                ML_tag = tuple_single[1].tolist()
        if MM_tag == '':
            return [], [], [], []
        MM_tag_split = MM_tag.split(";")
        hmC_tag = []
        hmC_percent = []
        mC_tag = []
        mC_percent = []
        mC_tag_head = ''
        hmC_tag_head = ''
        for MM_split in MM_tag_split:
            c_counting = MM_split.split(',')
            len_bases = len(c_counting[1:])
            percent = ML_tag[:len_bases]
            ML_tag = ML_tag[len_bases:]
            if 'm' in c_counting[0].lower() and 'c' in c_counting[0].lower():
                mC_tag_head = c_counting[0]
                mC_tag = c_counting[1:]
                mC_percent = percent[:]
            elif 'h' in c_counting[0].lower() and 'c' in c_counting[0].lower():
                hmC_tag_head = c_counting[0]
                hmC_tag = c_counting[1:]
                hmC_percent = percent[:]
        return mC_tag_head, hmC_tag_head, mC_tag, hmC_tag, mC_percent, hmC_percent

    def trimming_methyl_tuple(self, start, end):
        c_coors = []
        for i in range(self.sam_obj.query_length):
            if self.sam_obj.query_sequence[i] == 'C' or self.sam_obj.query_sequence[i] == 'c':
                c_coors.append(i)
        mC_tag_head, hmC_tag_head, mC_tag, hmC_tag, mC_percent, hmC_percent = self.unstring_methyl_tuple()
        mC_tag_coor = []
        mC_coor_count = 0
        for mC in mC_tag:
            mC_coor_count = mC_coor_count + int(mC)
            mC_tag_coor.append(c_coors[mC_coor_count])
            mC_coor_count = mC_coor_count + 1
        hmC_tag_coor = []
        hmC_coor_count = 0
        for hmC in hmC_tag:
            hmC_coor_count = hmC_coor_count + int(hmC)
            hmC_tag_coor.append(c_coors[hmC_coor_count])
            hmC_coor_count = hmC_coor_count + 1

        mC_tag_coor_trim = mC_tag_coor[:]
        for mC in mC_tag_coor:
            if mC < start: 
                mC_tag_coor_trim.pop(0)
                mC_percent.pop(0)
            elif mC >= end:
                mC_tag_coor_trim.pop(-1)
                mC_percent.pop(-1)
        for index in range(len(mC_tag_coor_trim)):
            mC_tag_coor_trim[index] = mC_tag_coor_trim[index] - start

        hmC_tag_coor_trim = hmC_tag_coor[:]
        for hmC in hmC_tag_coor:
            if hmC < start: 
                hmC_tag_coor_trim.pop(0)
                hmC_percent.pop(0)
            elif hmC >= end:
                hmC_tag_coor_trim.pop(-1)
                hmC_percent.pop(-1)
        for index in range(len(hmC_tag_coor_trim)):
            hmC_tag_coor_trim[index] = hmC_tag_coor_trim[index] - start
        
        middle_sequence = self.sam_obj.query_sequence[start:end]
        middle_c = []
        for index, character in enumerate(middle_sequence):
            if character == 'c' or character == 'C':
                middle_c.append(index)
        mC_counter = 0
        hmC_counter = 0
        mC_index_counter = 0
        hmC_index_counter = 0
        mC_tag_trimmed = []
        hmC_tag_trimmed = []
        for c in middle_c:
            if mC_index_counter < len(mC_tag_coor_trim) and c == mC_tag_coor_trim[mC_index_counter]:
                mC_tag_trimmed.append(mC_counter)
                mC_index_counter = mC_index_counter + 1
                mC_counter = 0
            else:
                mC_counter = mC_counter + 1
            if hmC_index_counter < len(hmC_tag_coor_trim) and c == hmC_tag_coor_trim[hmC_index_counter]:
                hmC_tag_trimmed.append(hmC_counter)
                hmC_index_counter = hmC_index_counter + 1
                hmC_counter = 0
            else:
                hmC_counter = hmC_counter + 1
        self.string_methyl_tuple(mC_tag_head, hmC_tag_head, mC_tag_trimmed, hmC_tag_trimmed, mC_percent, hmC_percent)

    def string_methyl_tuple(self, mC_tag_head, hmC_tag_head, mC_tag_trimmed, hmC_tag_trimmed, mC_percent, hmC_percent):
        mC_tag_string = mC_tag_head + "," + ",".join(map(str, mC_tag_trimmed)) + ";"
        hmC_tag_string = hmC_tag_head + "," + ",".join(map(str, hmC_tag_trimmed)) + ";"
        MM_string = hmC_tag_string + mC_tag_string
        ML_array = hmC_percent + mC_percent
        if ML_array:
            self.sam_obj.set_tag("MM", MM_string)
            self.sam_obj.set_tag("ML", ML_array)

    def alignment_eraser(self):
        self.sam_obj.cigarstring = None
        self.sam_obj.flag = 4
        self.sam_obj.reference_id = -1
        self.sam_obj.reference_start = -1
        self.sam_obj.mapping_quality = 0
        self.sam_obj.template_length = len(self.sam_obj.query_sequence)



    