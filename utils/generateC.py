#!/usr/bin/python

#Load required modules
import sys, os, argparse, math, random, time


def parse_arguments():
	"""
	Parse command line arguments
	Returns:
	"""

	description = "Generate a matrix C of tumor profiles from a random phylogeny."
	parser = argparse.ArgumentParser(description=description)
        parser.add_argument("-c", "--chromosomes", required=False, default=1, type=int, help="The number of chromosomes")
        parser.add_argument("-k", "--leaves", required=True, type=int, help="The number of leaves composing the rows of C")
        parser.add_argument("-n", "--segments", nargs='+', required=True, type=int, help="The number of profile segments, i.e. the columns of C")
        parser.add_argument("-u", "--events", nargs='+', required=True, type=int, help="The maximum number of events on each edge")
        parser.add_argument("-s","--seed", required=False, type=int, default=time.time(), help="Random seed")
        parser.add_argument("-d","--delproportion", required=False, type=float, default=0.4, help="The proportion of deletions than amplifications")
        parser.add_argument("-o","--directory", required=False, type=str, default="./", help="The directory of the output file")
	args = parser.parse_args()

        num_chromosomes=args.chromosomes
	num_leaves=args.leaves
	num_segments=args.segments
        max_events=args.events
        seed=args.seed
        del_proportion=args.delproportion
        directory=args.directory

        if (not(num_chromosomes) or not(num_leaves > 0)):
            raise ValueError('number of chromosomes, leaves, segments, and events must be greater than zero')

        if (len(num_segments) == 1):
                num_segments = [num_segments[0] for i in range(num_chromosomes)]
        elif (len(num_segments) != num_chromosomes):
                raise ValueError('the number of segments must be one or specified for each chromosome')

        if (len(max_events) == 1):
                max_events = [max_events[0] for i in range(num_chromosomes)]
        elif (len(max_events) != num_chromosomes):
                raise ValueError('the maximum number of events must be one or specified for each chromosome')

        if (del_proportion < 0.0 or del_proportion > 1.0):
            raise ValueError('the proportions of deletions vs amplifications must be a value between 0.0 and 1.0')

	return num_chromosomes, num_leaves, num_segments, max_events, seed, del_proportion, directory




class Profile:
    """
    """
    def __init__(self):
        self.name = -1
        self.profile = []
        self.right = None
        self.left = None

    def __repr__(self):
            if(self.left != None and self.right!=None):
                    return '%d-%s:L(%s)R(%s)' % (self.name, self.profile, self.left, self.right)
            else:
                    return '%d-%s' % (self.name, self.profile)

    def set_children(self, left, right):
        self.left = left
        self.right = right




def mutate_profile(profile, k, del_proportion, chro, events, source, target):
        result = list(profile)
        num_segments = len(result)
        num_events = random.randint(1, k);
        for i in range(num_events):
                interval = []
                interval.append(random.randint(0, num_segments-1))
                interval.append(random.randint(0, num_segments-1))
                interval.sort();
                event = 1
                if(random.random() <= del_proportion):
                        event = -1
                events[(chro, source, target, interval[0], interval[1])] += event
                for i in range(interval[0], interval[1]+1):
                        if(result[i] > 0):
                                result[i] = result[i] + event
        return result



def build_trees(num_leaves):
    result = []
    if(num_leaves == 1):
        result.append(Profile())
    else:
        for split in range(1, num_leaves):
            left_trees = build_trees(split)
            right_trees = build_trees(num_leaves - split)

            for l in left_trees:
                for r in right_trees:
                    subtree = Profile()
                    subtree.set_children(l, r)

                    result.append(subtree)
    return result



def enumerate_nodes(tree, counter):
        tree.name = counter
        counter += 1
        counter = enumerate_internals(tree.left, counter)
        counter = enumerate_internals(tree.right, counter)
        counter = enumerate_leaves(tree.left, counter)
        counter = enumerate_leaves(tree.right, counter)



def enumerate_internals(tree, counter):
        if(tree.left != None):
                tree.name = counter
                counter += 1

                counter = enumerate_internals(tree.left, counter)
                counter = enumerate_internals(tree.right, counter)
        return counter



def enumerate_leaves(tree, counter):
        if(tree.left == None):
                tree.name = counter
                counter += 1
        else:
                counter = enumerate_leaves(tree.left, counter)
                counter = enumerate_leaves(tree.right, counter)
        return counter



def label_tree(tree, num_segments, max_events, del_proportion, chro, events):
        tree.profile = [2 for i in range(num_segments)]

        tree.left.profile = mutate_profile(tree.profile, max_events, del_proportion, chro, events, tree.name, tree.left.name)
        label_tree_rec(tree.left, num_segments, max_events, del_proportion, chro, events)

        tree.right.profile = mutate_profile(tree.profile, max_events, del_proportion, chro, events, tree.name, tree.right.name)
        label_tree_rec(tree.right, num_segments, max_events, del_proportion, chro, events)



def label_tree_rec(node, num_segments, max_events, del_proportion, chro, events):
        if(node.left != None):
                node.left.profile = mutate_profile(node.profile, max_events, del_proportion, chro, events, node.name, node.left.name)
                label_tree_rec(node.left, num_segments, max_events, del_proportion, chro, events)

                node.right.profile = mutate_profile(node.profile, max_events, del_proportion, chro, events, node.name, node.right.name)
                label_tree_rec(node.right, num_segments, max_events, del_proportion, chro, events)


def extract_vertices(tree, chro, internals, leaves):
        if(tree.left == None):
                leaves[(chro, tree.name)] = tree.profile
        else:
                internals[(chro, tree.name)] = tree.profile
                extract_vertices(tree.left, chro, internals, leaves)
                extract_vertices(tree.right, chro, internals, leaves)


def print_tree(tree, f):
        if(tree.left != None):
                f.write(str(tree.name) + " -> " + str(tree.left.name) + "\n")
                f.write(str(tree.name) + " -> " + str(tree.right.name) + "\n")
                print_tree(tree.left, f)
                print_tree(tree.right, f)


def print_events(f, events, num_chromosomes, num_leaves, num_segments):
        for chro in range(num_chromosomes):
                for i in range(num_leaves - 1):
                        for j in range(i + 1, 2*num_leaves - 1):
                                for s in range(num_segments[chro]):
                                        for t in range(s, num_segments[chro]):
                                                if(events[(chro, i, j, s, t)] != 0):
                                                        f.write(str(chro) + " " + str(i) + " " + str(j)  + " " + str(s)  + " " + str(t)  + " " + str(events[(chro, i, j, s, t)]) + "\n");




def print_full_output(tree, internals, leaves, num_chromosomes, num_leaves, num_segments, events, max_events, seed, del_proportion, directory):
        filename = "simC_c" + str(num_chromosomes) + "_k" + str(num_leaves) + "_n"
        for num in num_segments:
                filename = filename + str(num)
        filename = filename + "_u"
        for k in max_events:
                filename = filename + str(k)
        filename = filename + "_s" + str(seed) + "_del" + str(del_proportion).replace(".", "") + ".true"

	path = os.path.join(directory, filename)
        with open(path, 'w') as f:
                f.write("#PARAMS\n")
                f.write(str(num_chromosomes) + " #number of chromosomes\n")
                f.write(str(num_leaves) + " #number of leaves\n")
                f.write(' '.join([str(i) for i in num_segments]) + " #number of segments for each chromosome\n")
                # f.write(' '.join([str(i) for i in max_events]) + " #maximum number of edge events for each chromosome\n")
                # f.write(str(seed) + " #seed\n")

                f.write("#PROFILES\n")

                for inner in range(num_leaves - 1):
                        f.write(str(inner) + " : ")
                        for chro in range(num_chromosomes):
                                f.write(' '.join([str(i) for i in internals[(chro, inner)]]))
                                if(chro < (num_chromosomes - 1)):
                                        f.write(" | ")
                        f.write("\n")


                for leaf in range(num_leaves):
                        id_leaf = leaf + num_leaves - 1;
                        f.write(str(id_leaf) + " : ")
                        for chro in range(num_chromosomes):
                                f.write(' '.join([str(i) for i in leaves[(chro, id_leaf)]]))
                                if(chro < (num_chromosomes - 1)):
                                        f.write(" | ")
                        f.write("\n")

                f.write("#EDGES\n")
                print_tree(tree, f)
                f.write("#EVENTS\n")
                print_events(f, events, num_chromosomes, num_leaves, num_segments)
        print 'OUTPUT WRITTEN IN:  %s' % (path)



def print_profile_output(tree, leaves, num_chromosomes, num_leaves, num_segments, max_events, seed, del_proportion, directory):
        filename = "simC_c" + str(num_chromosomes) + "_k" + str(num_leaves) + "_n"
        for num in num_segments:
                filename = filename + str(num)
        filename = filename + "_u"
        for k in max_events:
                filename = filename + str(k)
        filename = filename + "_s" + str(seed) + "_del" + str(del_proportion).replace(".", "") + ".input"

	path = os.path.join(directory, filename)
        with open(path, 'w') as f:
                f.write("#PARAMS\n")
                f.write(str(num_chromosomes) + " #number of chromosomes\n")
                f.write(str(num_leaves) + " #number of leaves\n")
                f.write(' '.join([str(i) for i in num_segments]) + " #number of segments for each chromosome\n")
                # f.write(' '.join([str(i) for i in max_events]) + " #maximum number of edge events for each chromosome\n")
                # f.write(str(seed) + " #seed\n")

                f.write("#PROFILES\n")
                for leaf in range(num_leaves):
                        id_leaf = leaf + num_leaves - 1;
                        f.write(str(id_leaf) + " : ")
                        for chro in range(num_chromosomes):
                                f.write(' '.join([str(i) for i in leaves[(chro, id_leaf)]]))
                                if(chro < (num_chromosomes - 1)):
                                        f.write(" | ")
                        f.write("\n")

        print 'OUTPUT WRITTEN IN:  %s' % (path)



def init_events(events, chro, num_leaves, num_segments):
        for i in range(num_leaves - 1):
                for j in range(i + 1, 2*num_leaves - 1):
                        for s in range(num_segments):
                                for t in range(s, num_segments):
                                        events[(chro, i, j, s, t)] = 0


# Main method
def main():

	#Load arguments
	num_chromosomes, num_leaves, num_segments, max_events, seed, del_proportion, directory = parse_arguments()
        random.seed(seed)

        trees = build_trees(num_leaves)
        # for i in range(len(trees)):
        #     print '%i -- %s' % (i, trees[i])

        tree = trees[random.randint(0, len(trees)-1)]
        enumerate_nodes(tree, 0)
        internals = {}
        leaves = {}
        events = {}

        for chro in range(num_chromosomes):
                init_events(events, chro, num_leaves ,num_segments[chro])
                label_tree(tree, num_segments[chro], max_events[chro], del_proportion, chro, events)
                extract_vertices(tree, chro, internals, leaves)

        print_full_output(tree, internals, leaves, num_chromosomes, num_leaves, num_segments, events, max_events, seed, del_proportion, directory)
        print_profile_output(tree, leaves, num_chromosomes, num_leaves, num_segments, max_events, seed, del_proportion, directory)



if __name__ == '__main__':
	main()
