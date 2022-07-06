# -*- coding: utf-8 -*-

from threading import Thread
import re
from os import listdir, path
import numpy as np
from itertools import chain

# Abstrakty - to jest potrzebne bardziej do 
# rozwijania roznych sposobow przetwarzania danych 
# az do otrzymania struktury z ktorej da sie policzyc genus


class GenusData(object):
    def __init__(self, name, mini, maxi, genuses):
        self.name = str(name)
        self.mini = int(mini)
        self.maxi = int(maxi)
        self.genuses = genuses


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
    def generate_inputs(cls, config):
        ''' abstract of generate_inputs method
        All you need is config dictionary with keys:
        directory - path for your data directory
        end_file - type of your file
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
    def create_workers(cls, input_class, config):
        ''' workers creator'''
        workers = []
        for input_data in input_class.generate_inputs(config):
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
            raise Exception("You don't have cleaned data")

    @staticmethod
    def przelicznik(tab):
        ''' remove empty spots'''
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

    
# wielowatkowosc - ciekawe czy dziala ? 
class FactorizeThread(Thread):
    ''' multi thread for Worker.map() method'''
    def __init__(self, worker):
        super().__init__()
        self.worker = worker

    def run(self):
        self.worker = self.worker.map()

#======================================
# Two columns data analysis - Proteins
#======================================

# analizy dla protein 
# dane tylko w dwoch kolumnach


class ProteinData(InputData):
    '''Load files from directory and 
    made two columns structure. 
    Don't use this class separately. 

    Parameters
    ==========
    path_dir: string
    
    Return
    ======
    path: string
    name: string, file name from path
    '''

    def __init__(self, path_dir):
        super().__init__()
        self.path = path_dir
        self.name = self.__take_name_from_path(self.path)

    def read(self):
        ''' read two columns data by ordinary open method.
        data could be separated by tab, |, ',' and -.
        '''
        data = []
        result = []
        with open(self.path) as file:
            for _, line in enumerate(file):
                if line.strip():
                    line = re.sub('[;:\t|,-]', ' ', line)
                    tokens = line.split()
                    if len(tokens) == 2:
                        result = [int(x) for x in tokens]
                        data.append(result)
        return data

    @staticmethod
    def __take_name_from_path(path_dir):
        '''
        helper method. 
        Take name of file from path
        '''
        path_split = path_dir.split('/')
        return path_split[len(path_split) - 1].split('.')[0]

    @classmethod
    def generate_inputs(cls, config):
        '''read all files in dir with some ends'''
        data_dir = config['directory']
        end = config['end_file']
        start = config['start_file'] if 'start_file' in config else None
        if start is not None:
            for name in [f for f in listdir(data_dir) if f.startswith(start)]:
                yield cls(path.join(data_dir, name))
        else:
            for name in [f for f in listdir(data_dir) if f.endswith(end)]:
                yield cls(path.join(data_dir, name))

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['path_dir']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)


class ProteinDataException(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class ProteinEmptyDataException(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class ProteinDataChimeraFormat(ProteinData):

    def read(self):
        ''' read two columns data from a Chimera format output line i.e.:
        ALA 8.A N        THR 11.A O       no hydrogen  2.876  N/A
        '''
        data = []
        result = []
        try:
            with open(self.path) as file:
                header = True
                for _, line in enumerate(file):
                    if header:
                        end_of_header = line.startswith('H-bonds')
                        if end_of_header:
                            header = False
                        continue
                    if line.strip():
                        line = re.sub('[\t|,-]', ' ', line)
                        tokens = line.split()
                        if len(tokens) > 4:
                            try:
                                result = [filter_token(tokens[1]), filter_token(tokens[4])]
                                data.append(result)
                            except ProteinDataException:
                                pass
        except ValueError as ve:
            raise ProteinDataException("Problem reading Chimera file: " + self.path + "\n" + str(ve))
        if len(data) == 0:
            raise ProteinEmptyDataException("Empty data for : " + self.path)
        return data


def filter_token(token):
    tok = token.split(".")[0]
    if tok.find("A") > 0:
        tok = tok.replace("A", "")
    try:
        return int(tok)
    except ValueError:
        raise ProteinDataException("Undesired contact")


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
        chords = len(array)
        if chords == 0:
            raise ProteinEmptyDataException("Protein List empty after cleaning")
        maxi = int(array.max())
        mini = int(array.min())
        length = maxi - mini + 1
        return mini, maxi, length, chords

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

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)


# dla analiz generujesz już tylko Workera
# z określeniem dla niego metody map 
class ProteinWorker(Worker):
    ''' class for proteine genus analysis'''

    def map(self):
        ''' map work for one protein '''
        self.name = self.input_data.name
        self.cl_data = ProteinList(self.input_data.read()).clean()
        self.mini, self.maxi,\
        self.length, self.nr_chord = self.cl_data.compute_info()
        self.genus = self.compute_genus()

    def info(self):
        return self.__dict__


# dla roznych analiz dodajesz tylko nową klasę z Workerem 
# i metodą map
class ProteinAnalysisWorker(Worker):
    ''' class for devided proteine analysis'''

    def map(self):
        ''' map work for one protein '''
        self.name = self.input_data.name
        self.cl_data = ProteinList(self.input_data.read()).clean()
        self.mini, self.maxi, \
        self.length, self.nr_chord = self.cl_data.compute_info()
        self.no_biff = ProteinList(self.remove_biff(self.cl_data)).change()
        self.b1_b2_list = self.__multiplicity(self.cl_data)
        self.genuses = self.domains(self.divide_and_compute())

    def __multiplicity(self, tab):
        '''compute biffurcation numbers'''
        b1_b2_list = []
        for i in range(self.mini, self.maxi + 1):
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

    def main_trace(self) -> list:
        return [self.mini, self.one_trace(self.no_biff)]

    def one_trace(self, tab: ProteinList, mini_index: int = 1) -> list:
        data_len = 0
        genuses = []
        b1_b2_sum = np.array(self.b1_b2_list[:mini_index]).sum()
        mini = np.array(self.b1_b2_list[:mini_index]).sum()
        for ix, el in enumerate(self.b1_b2_list[mini_index:]):
            data = []

            if not el:
                # liczba anionów w adresie = 0
                if not ix:
                    # jeśli jesteś na pierwszym elemencie
                    genuses.append(0)
                else:
                    genuses.append(genuses[-1])
            else:  # są jakieś aniony licz
                b1_b2_sum += el
                # wygenerowanie listy wiazan
                data = [x for x in tab if x[1] <= b1_b2_sum and \
                        x[1] >= mini and \
                        x[0] <= b1_b2_sum \
                        and x[0] >= mini]
                if data:
                    if data_len != len(data):
                        data_len = len(data)
                        data = self.przelicznik(data)
                        genuses.append(self.genus_one_backbone(data))
                    else:
                        genuses.append(genuses[-1])
                else:
                    genuses.append(0)
        return genuses

    def divide_and_compute(self):
        '''Divide data and compute genus for all data with lag 1'''
        matrix = []
        for ix, i in enumerate(range(self.mini, self.maxi)):
            matrix.append([i, self.one_trace(self.no_biff, ix + 1)])
        return matrix

    @staticmethod
    def matrix_plot(tab, mini=0):
        """convert genuses list to list for 3d plot"""
        result = []
        matrix = np.zeros(shape=(len(tab), len(tab)))
        for el in tab:

            for i, e in enumerate(el[1]):
                temp = []
                temp.append(el[0])
                temp.append(el[0] + i + 1)
                temp.append(e)
                matrix[el[0] - mini + i, el[0] - mini] = e
                result.append(temp)
        return matrix

    @staticmethod
    def __short_list(l: list = []) -> list:
        result = []
        s = max(l)
        for el in range(s + 1):
            result.append(l.count(el))
        return result

    def domains(self, l: list = []) -> list:
        result = []
        for el in l:
            temp = []
            temp.append(el[0])
            temp.append(self.__short_list(el[1]))
            result.append(temp)
        return result

    @staticmethod
    def genus_list(l: list = []) -> list:
        result = []
        for el in l:
            temp = []
            g = 0
            for i in el[1]:
                for j in range(i):
                    temp.append(g)
                g = g+1
            result.append([el[0], temp])
        return result

    def total_genus(self):
        return len(self.genuses[0][-1])-1

    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        try:
            del state['input_data']
        except:
            pass
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

# to na koncu


def data_analysis(worker_class, input_class, config):
    ''' helper for Thread with map() method '''
    threads = []
    workers = worker_class.create_workers(input_class, config)

    if config['multithreading']:
        for worker in workers:
            thread = FactorizeThread(worker)
            thread.start()
            threads.append(thread)
        for thread in threads:
            thread.join()
    else:
        for worker in workers:
            worker.map()
    return workers        


def get_input_class(config):
    if config['end_file'] == 'chim':
        return ProteinDataChimeraFormat
    else:
        return ProteinData


def protein_genus(config):
    '''run proteins analysis from directory'''
    return data_analysis(ProteinWorker, get_input_class(config), config)

'''
    Deviding two columns and compute genus from min_length
    to max_length by one lag. Return object
'''
def protein_genus_trace(config):
    '''run proteins deviding and analysis'''
    return data_analysis(ProteinAnalysisWorker, get_input_class(config), config)