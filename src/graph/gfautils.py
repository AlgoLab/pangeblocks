import gfapy

class PrintListGFA:
    def __init__(self):
        self.objs = []

    def append(self, obj):
        if type(obj) == gfapy.line.edge.link.link.Link:
            self.objs.append((obj.from_name, obj.to_name))
        else:
            self.objs.append(obj.name)
    
    def __contains__(self, obj):
        if type(obj) == gfapy.line.edge.link.link.Link:
            return (obj.from_name, obj.to_name) in self.objs
        else:
            return obj.name in self.objs

class SetGFA():
    def __init__(self):
        self.ix_ss = set()
        self.ix_ls = set()
        self.segments = dict()
        self.links = dict()
    
    def add(self, obj):
        if type(obj) == gfapy.line.edge.link.link.Link:
            k = (obj.from_name, obj.to_name)
            if not k in self.ix_ls:
                self.ix_ls.add(k)
                self.links[k] = obj
        else:
            k = obj.name
            if not k in self.ix_ss:
                self.ix_ss.add(k)
                self.segments[k] = obj

    def sort(self):
        for k in sorted(self.ix_ss):
            yield self.segments[k]

        for k in sorted(self.ix_ls):
            yield self.links[k]

def _rec_gn(gfa, s_name, radius, printset, full, ignore_back, ignore_forward):
    if radius == 0:
        seg = gfa.segment(s_name)
        printset.add(seg)
        return
    else:
        seg = gfa.segment(s_name)
        printset.add(seg)
        _fl = list()
        _bl = list()
        for l in gfa.dovetails:
            if l.from_segment == seg:
                _fl.append(l)
            elif l.to_segment == seg:
                _bl.append(l)
        
        igns = [False, False]
        if not ignore_back:
            for l in _bl:
                printset.add(l)
                if not full:
                    igns = [False, True]
                _rec_gn(gfa, gfa.segment(l.from_name), radius-1, printset, full, *igns)
        if not ignore_forward:
            for l in _fl:
                printset.add(l)
                if not full:
                    igns = [True, False]
                _rec_gn(gfa, gfa.segment(l.to_name), radius-1, printset, full, *igns)


def get_neighbourhood(gfa, s_name, radius=1, full=False):
    pset = SetGFA()
    _rec_gn(gfa, s_name, radius, pset, full, False, False)

    for p in pset.sort():
        print(p)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='GFA utils')
    subparsers = parser.add_subparsers(help='sub-command help', dest="subparser_name")

    sparser_nh = subparsers.add_parser('nh', help='Get the neighbordhood centered on NAME with RADIUS')
    sparser_nh.add_argument('--gfa', action='store', required=True, help='GFA file')
    sparser_nh.add_argument('--name', action='store', required=True, type=str, help='Name of the centered segment')
    sparser_nh.add_argument('--radius', action='store', default=1, type=int, help='Radius of the neighborhood [Default=1]')
    sparser_nh.add_argument('--full', action='store_true', help='Full neighborhood [Default=False]')

    args = parser.parse_args()

    if args.subparser_name == 'nh':
        gfa = gfapy.Gfa.from_file(args.gfa)
        get_neighbourhood(gfa, args.name, radius=args.radius, full=args.full)