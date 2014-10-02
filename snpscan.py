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


class SnpInterval(object):

    def __init__(self, interval, size, slide):
        self.intrval_snps = {}
        self.chrom, coord = interval.split(':')
        self.start, self.stop = map(int, coord.split('-'))
        self.size = size
        self.slide = slide

    def get_num_intervals(self):
        return len(self.intrval_snps)

    def _create_interval(self):
        istart = self.start  # interval start
        istop = self.start + self.size  # interval stop
        while istop < self.stop and istart < self.stop:
            intrval = map(str, [self.chrom, istart, istop])
            intrval = pybedtools.create_interval_from_list(intrval)
            yield intrval
            istart += self.slide
            istop += self.slide

        intrval = map(str, [self.chrom, istart, self.stop])
        intrval = pybedtools.create_interval_from_list(intrval)
        yield intrval

    def scan(self, snpfile):
        for intrval in self._create_interval():
            i = '%s:%d-%d' % (intrval.chrom, intrval.start, intrval.stop)
            self.intrval_snps[i] = []
            for snp in snpfile.filter(self._is_in, intrval):
                self.intrval_snps[i].append(SnpRec(snp.fields))

    def _is_in(self, snp, intrval):
        if (snp.start >= intrval.start and
                snp.end <= intrval.end):
            return True
        else:
            return False

    def display_intervals(self):
        for i in self.intrval_snps.keys():
            print '%s, total SNPs=%d' % \
                    (i, len(self.intrval_snps[i]))


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

    snpint = SnpInterval(interval, 50, 100)
    snpint.scan(snpfile)
    snpint.display_intervals()


if __name__=='__main__':
    main()
