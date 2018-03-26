# temporary container


class DASequenceContainer():
    def __init__(self, start_id, n_dots, dot_array_definitions, sequence_methods):
        """len(methods) and defined the number of sequences
        simultaniously produced """

        if len(n_dots) != len(dot_array_definitions) or len(dot_array_definitions) != len(sequence_methods):
            raise RuntimeError("n_dots, dot_array_definitions, sequence_methods must be arrays of the same length")

        self._current_id = start_id
        self.dot_array_definitions = dot_array_definitions
        self.sequence_methods = sequence_methods
        self.n_dots = n_dots

        self._processes = [None] * len(sequence_methods)
        for c in range(len(sequence_methods)):
            self.start_making_sequence(c)

    def start_making_sequence(self, num):
        da = DotArray(n_dots=self.n_dots[num],
                              dot_array_definition=self.dot_array_definitions[num],
                              array_id=self._current_id)
        self._current_id += 1
        self._processes[num] = MakeDASequenceProcess(max_dot_array=da,
                                                     method=self.sequence_methods[num],
                                                     n_trails=4)

    def sequence_available(self, num):
        return self._processes[num].sequence_available.is_set()

    def wait_sequence_available(self, num):
        return self._processes[num].sequence_available.wait()

    def wait_all_sequences_available(self):
        for x in self._processes:
            x.sequence_available.wait()

    def get_error(self, num):
        return self._processes[num].error

    def pop_sequence(self, num):
        """returns sequence and starts making new one"""
        rtn = self._processes[num].da_sequence
        self.start_making_sequence(num)
        return rtn



def OLD_make_multiple_dot_array_sequences(n_versions,
                                      methods,
                                      n_dots,
                                      stimulus_area_radius,
                                      dot_diameter_mean,
                                      dot_diameter_range=None,
                                      dot_diameter_std=None,
                                      dot_colour=None,
                                      minium_gap=1,
                                      sqeeze_factor=.9,
                                      subfolder = "stimuli",
                                      extension = EXTENSION,
                                      n_processes = cpu_count(),
                                      override_existing_arrays=False,
                                      save_erroneous_arrays=False):
    """!!!! # FIXME
    returns also feedback

    """

    try:
        os.mkdir(subfolder)
    except:
        pass

    # make processes
    processes = []
    for cnt in range(n_versions):
        for method in methods:
            filename = os.path.join(subfolder, method+"-"+str(cnt)+ extension)
            if override_existing_arrays or not os.path.isfile(filename):
                print("preparing "+filename)
                da = DotArray(n_dots=n_dots, stimulus_area_radius=stimulus_area_radius,
                              dot_diameter_mean=dot_diameter_mean, dot_diameter_range=dot_diameter_range,
                              dot_diameter_std=dot_diameter_std, dot_colour=dot_colour, minium_gap=minium_gap)

                if sqeeze_factor<1 and sqeeze_factor>0:
                    da.fit_convex_hull_area(convex_hull_area=da.convex_hull_area * sqeeze_factor)

                error_queue = Queue()
                sequence_queue = Queue()
                p = Process(target = make_da_sequence, args=(da, method, sequence_queue, error_queue))
                processes.append((p, error_queue, sequence_queue) )

    # run processes
    running=[]
    while len(processes)>0:
        # remove dead processes from running
        for p in running:
            if not p.is_alive():
                running.remove(p)

        # move process to running and start if possible
        if len(running) < n_processes:
            running.append(processes.pop(0))
            running[-1].start()
        else:
            sleep(1)

    # wait running processes are to terminate
    for p in running:
        p.join()

    realign_error_list = []
    while not error_queue.empty():
        realign_error_list.append(error_queue.get_nowait())

    return realign_error_list


