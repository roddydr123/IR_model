import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
from tqdm import tqdm


def update_grid(grid, grid_size, p):

    # array of random indices for one sweep.
    ijs = np.random.randint(0, grid_size, size=(grid_size**2, 2))

    # array of taskc for becoming inactive and infecting neighbour.
    random_probs = np.random.rand(grid_size**2, 2)

    # do a single sweep
    for step in range(grid_size**2):

        # fetch indices of point to update
        i, j = ijs[step]

        # if active
        if grid[i, j] == 1:

            # become inactive with probability (1-p).
            if p < random_probs[step, 0]:
                grid[i,j] = 0

            # infect neighbour.
            else:
                # choose which neighbour randomly.
                if random_probs[step, 1] < 0.25:
                    grid[i - 1, j] = 1
                elif random_probs[step, 1] > 0.25 and random_probs[step, 1] < 0.5:
                    grid[(i + 1) % grid_size, j] = 1
                elif random_probs[step, 1] > 0.5 and random_probs[step, 1] < 0.75:
                    grid[i, (j + 1) % grid_size] = 1
                else:
                    grid[i, j - 1] = 1

    return grid


def animation(grid, grid_size, p, vis):

    # if we wanted a visualisation, make a figure for it.
    if vis:
        fig, ax = plt.subplots()
        im = ax.imshow(grid, animated=True, cmap="gray")
        # cbar = fig.colorbar(im, ax=ax)

    infected_list = []
    nsteps = 10000

    for i in range(nsteps):

        # move one step forward in the simulation, updating at every point.
        grid = update_grid(grid, grid_size, p)

        infected_list.append(np.sum(grid == 1)/grid_size**2)

        # every 50 sweeps update the animation.
        if i % 10 == 0 and vis:
            
            plt.cla()
            im = ax.imshow(grid, interpolation=None, animated=True, cmap="gray")
            plt.draw()
            plt.pause(0.00001)

    if not vis:
        plt.plot(range(nsteps), infected_list)
        plt.ylim(0, 1)
        plt.show()


def taskb(grid, grid_size, p):

    infected_list = []
    nsteps = 200

    for i in tqdm(range(nsteps)):

        # move one step forward in the simulation, updating at every point.
        grid = update_grid(grid, grid_size, p)

        infected_list.append(np.sum(grid == 1)/grid_size**2)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(range(nsteps), infected_list)
    ax1.set_title(f"fraction of alive cells over time for p = {p}")
    ax1.set_ylabel("A/N")
    ax1.set_xlabel("sweeps")
    ax1.set_ylim(0, 1)
    plt.show()


def taskc():
    """Compute the average fraction of active sites."""
    p_list = np.arange(0.55, 0.7, 0.005)
    grid_size = 50

    nsteps = 1000

    average_infected_list = []

    for p in tqdm(p_list):
        grid = np.random.randint(2, size=(grid_size, grid_size))
        infected_list = []

        for i in range(nsteps):
            grid = update_grid(grid, grid_size, p)

            if i > 100:
                total_infected = np.sum(grid == 1)
                if total_infected != 0:
                    infected_list.append(total_infected/grid_size**2)
                else:
                    infected_list += [0] * (nsteps - i)
                    break

        average_infected_list.append(np.average(infected_list))

    plt.plot(p_list, average_infected_list)
    plt.title(f"average survival in {nsteps} steps versus p")
    plt.ylabel("average # infected / N")
    plt.xlabel("sweeps")
    plt.show()


def get_jacknife_error(data, interesting_quantity_func):
    """
    Working jacknife error function.
    
    data: raw data from which interesting quantity is calculated.
    
    interesting quantity func: function which takes raw data and returns interesting quantity.
    """

    # make sure data to be resampled is an array.
    data = np.array(data)

    # calculate the interested quantity using all data.
    all_samples = interesting_quantity_func(data)

    # prepare an array for the quantities calculated with resampling.
    resampled_quantity = np.zeros(len(data))

    # loop over all the data and remove one sample each time.
    for i in range(len(data)):

        # array with all but one data point in it.
        resample = data[np.arange(len(data)) != i]

        # calculate the interesting quantity with the slightly reduced dataset.
        resampled_quantity[i] = interesting_quantity_func(resample)

    # find the error on the interesting quantity using the calculated values.
    error = np.sqrt(np.sum((resampled_quantity - all_samples) ** 2))
    return error


def get_bootstrap_error(data, interesting_quantity_func):
        """working bootstrap method."""

        # how many resamples to do. 1000 should work well.
        k = 1000

        # prepare array for quantities calculated with resampling.
        resampled_quantity = np.zeros(k)

        # resample k times.
        for i in range(k):

            # take a sample from the data, same length as the data set but
            # resampled with replacement.
            resample = np.random.choice(data, len(data))

            # calculate the interesting quantity with resampled dataset.
            resampled_quantity[i] = interesting_quantity_func(resample)

        # find the standard deviation of the newly calculated values.
        error = np.sqrt(np.var(resampled_quantity))
        return error


def taskd():
    p_list = np.arange(0.55, 0.7, 0.005)
    grid_size = 50

    nsteps = 1000

    variance_list = []
    errors = []

    for p in tqdm(p_list):
        grid = np.random.randint(2, size=(grid_size, grid_size))
        infected_list = []

        for i in range(nsteps):
            grid = update_grid(grid, grid_size, p)

            if i > 500:
                total_infected = np.sum(grid == 1)
                if total_infected != 0:
                    infected_list.append(total_infected)
                else:
                    infected_list += [0] * (nsteps - i)
                    break

        variance_list.append(np.var(infected_list) / grid_size**2)
        errors.append(get_jacknife_error(infected_list, np.var) / grid_size**2)

    plt.errorbar(p_list, variance_list, yerr=errors)
    plt.title(f"variance versus probability")
    plt.ylabel("variance")
    plt.xlabel("p")
    plt.show()


def taske():
    # survival probability over time.

    p_list = [0.6, 0.625, 0.65]
    grid_size = 50
    nsteps = 300
    repeats = 200

    overall_infected_over_time = []

    for p in p_list:

        ijs = np.random.randint(grid_size, size=(repeats, 2))
        p_survival = []

        # repeat with each p a few times.
        for repeat in tqdm(range(repeats)):

            # set a random cell to one.
            grid = np.zeros((grid_size, grid_size))
            i, j = ijs[repeat]
            grid[i,j] = 1

            infected_over_time = []

            # run one simulation.
            for n in range(nsteps):
                grid = update_grid(grid, grid_size, p)

                infected = np.sum(grid == 1)

                # stop the experiment if all cells die.
                if infected > 0:
                    infected_over_time.append(1)
                else:
                    infected_over_time += [0] * (nsteps - n)
                    break

            p_survival.append(infected_over_time)

        overall_infected_over_time.append(np.sum(p_survival, axis=0))

    for ydata, prob in zip(overall_infected_over_time, p_list):
        np.savetxt(f"{prob}inf1.txt", ydata)


def analyse_taske():
    p_list = [0.6, 0.625, 0.65]
    nsteps = 300
    for prob in p_list:
        ydata = np.loadtxt(f"{prob}inf1.txt")
        plt.plot(np.log10(range(nsteps)), np.log10(ydata), label=prob)
    plt.xlabel("log10(sweeps)")
    plt.ylabel("log10(number of simulations with surviving cells)")
    plt.legend()
    plt.show()


def main():
    """Evaluate command line args to choose a function.
    """


    mode = sys.argv[1]

    grid_size = 50

    if mode == "vis":
        p = float(sys.argv[2])
        grid = np.random.randint(2, size=(grid_size, grid_size))
        animation(grid, grid_size, p, True)
    elif mode == "b":
        p = float(sys.argv[2])
        grid = np.random.randint(2, size=(grid_size, grid_size))
        taskb(grid, grid_size, p)
    elif mode == "c":
        taskc()
    elif mode == "d":
        taskd()
    elif mode == "e":
        taske()
    elif mode == "ae":
        analyse_taske()
    else:
        print("wrong args")


main()
