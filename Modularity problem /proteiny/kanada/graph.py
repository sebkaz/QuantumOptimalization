from threading import Thread
import numpy as np


class ProteinList(list):
    '''Jedna z najlepszych klas jakie napisałem. 
    Bazuje na rozszerzeniu możliwości listy. 
    Cleaning operation for Your two columns data.
    '''
    def __init__(self, members):
        super().__init__(members)

    def numpy_array(self):
        return np.array(self)

    def compute_info(self):
        array = self.numpy_array()
        #chords = len(array)
        maxi = array.max()
        mini = array.min()
        length = maxi - mini + 1
        return array.min(), array.max(), length, len(array)

    def change(self):
        result = [[row[1], row[0]] if row[0] > row[
            1] else [row[0], row[1]] for row in self]
        return ProteinList(result)

    def rem_same_chords(self):
        result = [list(x) for x in set(tuple(x) for x in self)]
        return ProteinList(sorted(result))

    def rem_inout(self):
        return ProteinList([row for row in self if row[0] != row[1]])

    def zero_rem(self):
        result = self.numpy_array()
        if 0 in result:
            result = result + 1
        return ProteinList(result.tolist())

    def clean(self):
        return ProteinList(self.change().rem_same_chords().rem_inout().zero_rem())

#======================================
# Two columns data analysis - Graph
#======================================

class InputData(object):
    '''Abstract class of Input Data.
    You can generate different way to read 
    your data by read method, and take many 
    different types files with many structures.
    '''

    def read(self):
        '''Abstract of read method '''
        raise NotImplementedError

    @classmethod
    def generate_inputs(cls):
        ''' abstract of generate_inputs method
        '''
        raise NotImplementedError


class Worker(object):
    ''' Worker abstract class '''

    def __init__(self, input_data):
        self.input_data = input_data

    def map(self):
        ''' abstract map method'''
        raise NotImplementedError

    @classmethod
    def create_workers(cls, input_class, graphs, names):
        ''' workers creator'''
        workers = []
        for input_data in input_class.generate_inputs(graphs, names):
            workers.append(cls(input_data))
        return workers

    def compute_genus(self):
        ''' methods for genus computation from self.data
            You can use them only if You have self.cl_data
            Returns
            -------
            genus_ : int
                genus value for data

            nobiff_data_ : list with [in,out]
                data without biffurcations
        '''
        if hasattr(self, 'cl_data'):
            return self.genus_one_backbone(self.remove_biff(self.cl_data))
        else:
            raise "You don't have cleaned data yet"

    @staticmethod
    def przelicznik(tab):
        ''' remove empty spots'''
        from itertools import chain
        jeden = list(chain.from_iterable((row[0], row[1]) for row in tab))
        for i in range(1, max(jeden)):
            if i not in jeden and i < max(jeden):
                while i not in jeden:
                    for index, item in enumerate(jeden):
                        if item > i:
                            jeden[index] -= 1

        result = [[jeden[i], jeden[i + 1]] for i in range(0, len(jeden), 2)]
        return result

    def remove_biff(self, tab):
        return self.__resolve_bifurcations(self.przelicznik(tab), len(tab))
    
    @staticmethod
    def __bifurcations_part(tab, b_par, bif):
        ''' part of biff function '''
        while True:
            log = False
            for k in range(b_par - 1):
                if tab[bif[k]][1] < tab[bif[k + 1]][1]:
                    log = True
                    tab[bif[k]][1], tab[
                        bif[k + 1]][1] = tab[bif[k + 1]][1], tab[bif[k]][1]
            if not log:
                break
        return tab

    def __resolve_bifurcations(self, tab, n_chords):
        ''' Methods for resolve biffurcations in structure '''
        for i in range(1, 2 * n_chords):
            bif1 = []
            bif2 = []
            for j in range(n_chords):
                if tab[j][0] == i or tab[j][1] == i:
                    if tab[j][1] < i:
                        bif1.append(j)
                    elif tab[j][0] < i:
                        bif1.append(j)
                        tab[j][1] = tab[j][0]
                        tab[j][0] = i
                    else:
                        bif2.append(j)
            b1_len, b2_len = len(bif1), len(bif2)
            if b1_len > 1:
                tab = self.__bifurcations_part(tab, b1_len, bif1)
            if b2_len > 1:
                tab = self.__bifurcations_part(tab, b2_len, bif2)
            if b1_len + b2_len > 1:
                for row2 in tab:
                    if row2[0] > i:
                        row2[0] += b1_len + b2_len - 1
                    if row2[1] > i:
                        row2[1] += b1_len + b2_len - 1
            for k in range(1, b1_len):
                tab[bif1[k]][0] += k
            for l_n in range(b2_len):
                tab[bif2[l_n]][0] += b1_len + l_n
        return tab    
    @staticmethod
    def genus_one_backbone(tab):
        '''compute genus from clean data '''
        n_chords = len(tab)
        n_spots = 4 * n_chords
        remaining_chain = [1] * n_spots
        remaining_spots = n_spots
        edge = 0  # to chcemy policzyć - ilość brzegów
        while remaining_spots > 0:
            edge += 1
            nr_j = 0
            nr_k = 1
            while remaining_chain[nr_j] == 0:
                nr_j += 1
                if nr_j % 2 == 0:
                    nr_k += 1  # number of the atom
            side = nr_j % 2
            # print "Number: %i, Atom: %i, Side: %i " % (j, k, side)
            k_boundrary = nr_k
            side_boundrary = side
            remaining_chain[(k_boundrary - 1) * 2 + side_boundrary] = 0

            while True:
                i = 0
                while tab[i][0] != k_boundrary and tab[i][1] != k_boundrary:
                    i += 1
                if tab[i][0] == k_boundrary:
                    k_boundrary = tab[i][1]
                else:
                    k_boundrary = tab[i][0]
                side_boundrary = (1 + side_boundrary) % 2
                remaining_chain[(k_boundrary - 1) * 2 + side_boundrary] = 0
                if side_boundrary == 0:
                    k_boundrary = k_boundrary - 1
                else:
                    k_boundrary = k_boundrary + 1
                if k_boundrary == 0:
                    k_boundrary = 2 * n_chords
                if k_boundrary == 2 * n_chords + 1:
                    k_boundrary = 1
                side_boundrary = (1 + side_boundrary) % 2
                remaining_chain[(k_boundrary - 1) * 2 + side_boundrary] = 0
                remaining_spots = remaining_spots - 2

                if 2 * (k_boundrary - 1) + side_boundrary == 2 * (nr_k - 1) + side:
                    break

        genus = (2 + n_chords - 1 - edge) / 2
        # print "Number of boundaries: %i, genus: %i" % (r,genus)
        return int(genus)


# analizy dla grafów 
# dane tylko w dwoch kolumnach
# =================================
# Graph DATA
# =================================
class GraphData(InputData):
    '''Load list for grapf'''

    def __init__(self, graph: list, name: str):
        super().__init__()
        self.name = name
        self.data = graph


    def read(self):
        ''' read two columns list from graph generator
        '''
        return self.data
    
    @classmethod
    def generate_inputs(cls, graphs: list, names: list):
        '''read all files in dir with some ends'''
        for graph, name in zip(graphs, names):
            yield cls(graph, name)

# ==========================
# GRAPH WORKER
# ==========================
class GraphWorker(Worker):
    ''' class for proteine genus analysis'''

    def map(self):
        ''' map work for one protein '''
        self.name = self.input_data.name
        self.cl_data = ProteinList(self.input_data.read()).clean()
        self.mini, self.maxi,\
        self.length, self.nr_chord = self.cl_data.compute_info()
        self.genus = self.compute_genus()

# ============================
# GRAPH ANALYSIS WORKER
# ============================
class GraphAnalysisWorker(Worker):
    ''' class for devided graph analysis'''

    def map(self):
        ''' map work for one graph '''
        self.name = self.input_data.name
        self.cl_data = ProteinList(self.input_data.read()).clean()
        self.mini, self.maxi, \
        self.length, self.nr_chord = self.cl_data.compute_info()
        self.no_biff = ProteinList(self.remove_biff(self.cl_data)).change()
        self.b1_b2_list = self.__multiplicity(
            self.cl_data, self.mini, self.maxi)
        self.genuses = self.devide_and_compute(
            self.mini, self.maxi, self.no_biff)


    @staticmethod
    def __multiplicity(tab, minimum, maximum):
        '''compute biffurcation numbers'''
        b1_b2_list = []
        for i in range(minimum, maximum + 1):
            b1 = 0
            b2 = 0
            for j in range(len(tab)):
                if tab[j][0] == i or tab[j][1] == i:
                    if tab[j][1] < i:
                        b1 += 1
                    elif tab[j][0] < i:
                        b1 += 1
                    else:
                        b2 += 1
            b1_b2_list.append(b1 + b2)
        return b1_b2_list

    def devide_and_compute(self, mini, maxi, tab):
        '''Devide data and compute genus for all data with lag 1'''
        b1_b2_sum = 0
        counter = 0
        genuses = []
        data = []
        length_data = 0
        for element in range(mini, maxi + 1):
            b1_b2_sum += self.b1_b2_list[counter]
            counter += 1
            data = [x for x in tab if x[1] <= b1_b2_sum and x[0] <= b1_b2_sum]
            if data != []:
                if len(data) != length_data:
                    length_data = len(data)
                    data = self.przelicznik(data)
                    genuses.append(self.genus_one_backbone(data))
                else:
                    genuses.append(genuses[(len(genuses)-1)])
            else:
                genuses.append(0)
        return genuses


class FactorizeThread(Thread):
    ''' multi thread for Worker.map() method'''
    def __init__(self, worker):
        super().__init__()
        self.worker = worker

    def run(self):
        self.worker = self.worker.map()


def data_analysis(worker_class, input_class, graphs, nazwy):
    '''helper for Thread with map() method'''
    threads = []
    workers = worker_class.create_workers(input_class, graphs, nazwy)
    for worker in workers:
        thread = FactorizeThread(worker)
        thread.start()
        threads.append(thread)
    for thread in threads:
        thread.join()
    return workers     

def graph_genus(graphs, names):
    '''run proteins analysis from directory'''
    return data_analysis(GraphWorker, GraphData, graphs, names)

def graph_genus_trace(graphs, names):
    '''run graph deviding and analysis'''
    return data_analysis(GraphAnalysisWorker, GraphData, graphs, names)