import os;
import sys;

def edit_distance(s1, s2):
    if (len(s1) < len(s2)):
        return edit_distance(s2, s1);
    if (len(s2) == 0):
        return len(s1);
    previous_row = range(len(s2) + 1);
    for i, c1 in enumerate(s1):
        current_row = [i + 1];
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1;
            deletions = current_row[j] + 1;
            substitutions = previous_row[j] + (c1 != c2);
            current_row.append(min(insertions, deletions, substitutions));
        previous_row = current_row;
    return previous_row[-1];

class Dataset:
    def __init__(self, line=''):
        if (line == ''):
            self.reads_path = '';
            self.type = '';         ### nanopore/pacbio/single/paired/mate
            self.frag_len = 0;     ### Length of the fragment for paired end reads, insert size for mate pair libraries.
            self.frag_stddev = 0;   ### Standard deviation of the above length.
            self.reads_path_a = '';
            self.reads_path_b = '';
        else:
            split_line = line.split(',');

            self.type = split_line[0];
            if (self.type == 'paired' or self.type == 'mate'):
                if (len(split_line) < 5):
                    sys.stderr.write('ERROR: Five arguments need to be specified: "reads_type,reads_path_a,reads_path_b,frag_len,frag_stddev"!\n');
                    return;
                if (os.path.basename(split_line[1]) == ''):
                    sys.stderr.write('ERROR: Reads file path not correctly specified! Exiting.\n');
                    exit(1);
                if (os.path.basename(split_line[2]) == ''):
                    sys.stderr.write('ERROR: Reads file path not correctly specified! Exiting.\n');
                    exit(1);

                self.reads_path_a = os.path.abspath(split_line[1]);
                self.reads_path_b = os.path.abspath(split_line[2]);
                self.frag_len = int(split_line[3]);
                self.frag_stddev = int(split_line[4]);

                if (edit_distance(self.reads_path_a, self.reads_path_b) > 1):
                    sys.stderr.write('ERROR: Paths to paired-end/mate-pair libraries should differ only in 1 character! Exiting.\n');
                    exit(1);
                # print 'self.reads_path_a = %s' % (self.reads_path_a);
                wildcard_path = list(self.reads_path_a);
                # print wildcard_path;
                # print '';
                d = zip(self.reads_path_a, self.reads_path_b);
                # print d;

                current_char = 0;
                for i,j in d:
                    if (i != j):
                        wildcard_path[current_char] = '?';
                    current_char += 1;
                self.reads_path = ''.join(wildcard_path);

                # elif (self.type == 'nanopore' or self.type == 'pacbio' or self.type == 'single'):
            else:
                if (len(split_line) < 2):
                    sys.stderr.write('ERROR: Two arguments need to be specified: "reads_type,reads_path"!\n');
                    return;
                if (os.path.basename(split_line[1]) == ''):
                    sys.stderr.write('ERROR: Reads file path not correctly specified! Exiting.\n');
                    exit(1);
                self.reads_path = os.path.abspath(split_line[1]);

            # else:
            #     sys.stderr.write('ERROR: Unknown type of reads specified as parameter: "%s"!\n' % (line));
            #     return;
