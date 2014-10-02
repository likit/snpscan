#!/usr/bin/env python
'''Description'''
import sys
import pybedtools


class Snp(object):
    def __init__(self, chrom, pos, ref, alt,
                            score, genotype, sample):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.score = score
        self.genotype = genotype
        self.sample = sample

    def __str__(self):
        return 'sample=%d:[%s:%s, %s, %s, %s, %d]' % \
                    (self.sample, self.chrom,
                            self.pos, self.ref,
                            self.alt, self.genotype,
                            self.score)

    def __repr__(self):
        return self.__str__()


class SnpRec(object):
    def __init__(self, fields):
        self.records = []
        self.fields = fields
        self.chrom, self.pos, _, self.ref, \
                self.alt, self.score = self.fields[:6]

        sample = 1
        for f in fields[9:]:
            genotype = f.split(':')[0]
            self.records.append(Snp(self.chrom, self.pos,
                                    self.ref, self.alt,
                                    float(self.score), genotype, sample))
            sample += 1

    def __str__(self):
        return '%s:%s [%s/%s]' % (self.chrom, self.pos,
                                    self.ref, self.alt)

    def __repr__(self):
        return self.__str__()


class SnpInterval(object):

    def __init__(self, interval, size, slide):
        self.window_snps = {}
        self.chrom, coord = interval.split(':')
        self.start, self.stop = map(int, coord.split('-'))
        self.size = size
        self.slide = slide

    def get_num_windows(self):
        return len(self.window_snps)

    def _create_window(self):
        istart = self.start  # interval start
        istop = self.start + self.size  # interval stop
        while istop < self.stop and istart < self.stop:
            window = map(str, [self.chrom, istart, istop])
            window = pybedtools.create_interval_from_list(window)
            yield window
            istart += self.slide
            istop += self.slide

        window = map(str, [self.chrom, istart, self.stop])
        window = pybedtools.create_interval_from_list(window)
        yield window

    def scan(self, snpfile):
        for intrval in self._create_window():
            i = '%s:%d-%d' % (intrval.chrom, intrval.start, intrval.stop)
            self.window_snps[i] = []
            for snp in snpfile.filter(self._is_in, intrval):
                self.window_snps[i].append(SnpRec(snp.fields))

    def _is_in(self, snp, window):
        if (snp.start >= window.start and
                snp.end <= window.end):
            return True
        else:
            return False

    def display_windows(self):
        for i in self.window_snps.keys():
            print '%s, total SNPs=%d' % \
                    (i, len(self.window_snps[i]))


def main():
    '''Main function'''

    try:
        snpfile = sys.argv[1]
        interval = sys.argv[2]
        slide = sys.argv[3]
    except:
        print >> sys.stderr, 'Cannot open the input file.'
        raise SystemExit

    snpfile = pybedtools.BedTool(snpfile)  # pybedtool object

    # create an interval and scan for snps
    intrval = SnpInterval(interval, 50, 100)
    intrval.scan(snpfile)

    # display all sliding windows
    intrval.display_windows()

    # list all snps in each sliding window
    for rec in intrval.window_snps['chrM:50-100']:
        print rec
        for snp in rec.records[:2]:
            print snp,
        print '...'


if __name__=='__main__':
    main()
