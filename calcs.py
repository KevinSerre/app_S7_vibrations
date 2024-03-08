import math
from dataclasses import dataclass
from typing import Optional


# density of pedestrians (Ped / m**2) corresponding to each traffic class as per CSA S7 Table 7.1:
S7_PED_DENSITIES = {"TC1": 0.1, "TC2": 0.2, "TC3": 0.5, "TC4": 1.0, "TC5": 1.5}

# acceleration limits (m / s**2) corresponding to each comfort level as per CSA S7 Table 7.2
COMFORT_LEVELS_VERTICAL = {"Maximum": 0.5 , "Mean": 1.0, "Minimum": 2.5, "Unacceptable": float('inf')}
COMFORT_LEVELS_LATERAL = {"Maximum": 0.1 , "Mean": 0.5, "Minimum": 0.8, "Unacceptable": float('inf')}

# critical frequencies range
CRITICAL_VERTICAL_RANGE_START = 1    # Hz
CRITICAL_VERTICAL_RANGE_END = 13     # Hz

CRITICAL_LATERAL_RANGE_START = 0.3  # Hz
CRITICAL_LATERAL_RANGE_END = 6.5    # Hz
# average pedestrian weight in Newton
G0 = 700
# average pedestrian mass in kg
G0_mass = G0 / 10

# corresponding factor for the resonating harmonics 1 to 5: 
VERTICAL_DYNAMIC_LOAD_FACTOR = [0.4, 0.1, 0.42, 0.041, 0.027]

# corresponding factor for the resonating harmonics 1 to 5: 
LATERAL_DYNAMIC_LOAD_FACTOR = [0.05, 0.01, 0.023, 0.0043, 0.011] 


@dataclass
class PedestrianBridge:
    """
    A dataclass to describe the data required to check pedestrian-induced vibrations 
    of a pedestrian bridge as described in CSA S7.
    
    Assumptions: 
        - All values are in SI units
        - Pedestrian bridge is simply-supported
    """
    L: float                # length of bridge
    B: float                # width of bridge accessible to pedestrians
    damping: float          # structural damping value
    mass: float             # mass per unit length of bridge (structural + non-structural components)
    E: float                # modulus of elasticity of material
    I_vert: float           # second moment of area in the vertical direction
    I_lat: float            # second moment of area in the lateral direction
    G: Optional[float] = None                
    J: Optional[float] = None                # 
    mmi: Optional[float] = None              # mass moment of inertia
    

    def bending_frequencies(
        self,
    ) -> list[tuple[int, float]]:
        """
        Returns a list containing the natural bending frequencies of the bridge in the vertical direction.
        The first item is the fondamental vertical bending frequency, i.e first harmonic.
        """
        # calculate for the first 10-th harmonics
        harmonics = range(1, 11, 1)
        frequencies = [beam_bending_frequency(n, self.L, self.E, self.I_vert, self.mass) for n in harmonics] 

        return list( zip(harmonics, frequencies))


    def lateral_frequencies(
        self,
    ) -> list[tuple[int, float]]:
        """
        Returns a list containing the natural bending frequencies of the bridge in the lateral direction.
        The first item is the fondamental lateral bending frequency, i.e first harmonic.
        """
        # calculate for the first 10-th harmonics
        harmonics = range(1, 11, 1)
        frequencies = [beam_bending_frequency(n, self.L, self.E, self.I_lat, self.mass) for n in harmonics] 

        return list( zip(harmonics, frequencies))


    def torsion_frequencies(
        self,
    ) -> list[tuple[int, float]]:
        """
        Returns a list containing the natural torsional frequencies of the bridge.
        The first item is the fondamental torsional frequency, i.e first harmonic.
        """
        try:
            # calculate for the first 10-th harmonics
            harmonics = range(1, 11, 1)
            frequencies = [beam_torsion_frequency(n, self.L, self.G, self.J, self.mmi) for n in harmonics] 
            return list( zip(harmonics, frequencies))
        except:
            return []

    def calculate_max_vertical_acceleration(
        self,
        d: float,
        mass: float,
        n: int,
        f_n: float,
    ) -> float:
        """
        Returns the maximum vertical acceleration of the pedestrian bridge as per CSA S7 Cl.7.3.3.7.
        Assumes the bridge is a simply-supported beam with a distributed mass and a distributed harmonic loading along the length.
        """

        f_s = expected_vertical_walking_frequency(d)
        m = vertical_resonating_harmonic(f_n, f_s)
        n_eq = equivalent_number_of_pedestrians(d, self.L, self.B, self.damping)
        load = vertical_harmonic_load(m, n_eq) * self.B
        return beam_bending_max_acceleration(n, mass, self.damping, load)

    def calculate_max_lateral_acceleration(
        self,
        d: float,
        mass: float,
        n: int,
        f_n: float,
    ) -> float:
        """
        Returns the maximum lateral acceleration of the pedestrian bridge as per CSA S7 Cl.7.3.3.7.
        Assumes the bridge is a simply-supported beam with a distributed mass and a distributed harmonic loading along the length.
        """

        f_s = expected_lateral_walking_frequency(d)
        m = lateral_resonating_harmonic(f_n, f_s)
        n_eq = equivalent_number_of_pedestrians(d, self.L, self.B, self.damping)
        load = lateral_harmonic_load(m, n_eq) * self.B
        return beam_bending_max_acceleration(n, mass, self.damping, load)


# general functions : 

def beam_bending_frequency(
    n: int,
    L: float,
    E: float,
    I: float,
    m: float,
) -> float:
    """
    Returns the bending frequency of the n-harmonic for a simply-supported beam.

    'n': Rank of harmonic
    'L': Length of beam
    'E': Modulus of elasticity
    'I': Second moment of area in the direction of bending
    'm': Mass per unit length of the beam

    This function assumes consistent units are used for each parameter.
    """
    return n**2 * math.pi / (2 * L**2) * math.sqrt(E * I / m)


def beam_torsion_frequency(
    n: int,
    L: float,
    G: float,
    J: float,
    mmi: float,
) -> float:
    """
    Returns the bending frequency of the n-harmonic for a simply-supported beam.

    'n': Rank of harmonic
    'L': Length of beam
    'G': Shear modulus of elasticity
    'J': Torsional constant
    'mmi': Mass moment of inertia of the beam

    This function assumes consistent units are used for each parameter.
    """
    return n / (2 * L) * math.sqrt(G * J / mmi)


def beam_bending_max_acceleration(
    n: int,
    m: float,
    damping: float,
    F: float,
) -> float:
    """
    Returns the maximum acceleration of a simply-supported beam with distributed mass and subjected 
    to a distributed harmonic load having the same frequency as the n-th harmonic bending mode of the beam.

    'n': Rank of harmonic
    'm': Mass per unit length of the beam
    'damping': Structural damping value
    'F': Amplitude of the harmonic load per unit length of beam

    This function assumes consistent units are used for each parameter.
    """
    return 2 * F / (n * math.pi * m * damping)


# helper / utils functions :

def round_to_nearest_odd(
    n: float,
) -> int:
    """
    Returns an integer that is the nearest odd of 'n'.
    """
    lower_int = math.trunc(n)

    if lower_int % 2 == 1:
        nearest_odd = lower_int
    else:
        nearest_odd = math.trunc(n + 1)

    return nearest_odd


def reorder_dict(
    d: dict[str: dict[ str: list[ tuple[int, float, float]]]],
) -> dict[str: dict[ str: list[ tuple[int, float, float]]]]:
    """
    Returns a reordered nested dictionnary where second-level keys become first-level keys.
    Reordered first-level keys are the traffic classes defined in CSA S7.
    """

    new_dict = {key: {} for key in S7_PED_DENSITIES.keys()}
    for outer_key, inner_dict in d.items():
        for inner_key, value in inner_dict.items():
            new_dict[inner_key][outer_key] = value
    return new_dict


def merge_data(
    d: dict[str: dict[ str: list[ tuple[int, float, float]]]],
) -> dict[ str: list[ tuple[int, float, float]]]:
    """
    Returns a dictionnary where the values corresponding to second-level keys are merged together and sorted based on the 2nd item.
    2nd item of the tuples are the frequencies.
    """
    merged_data = {}
    for key, value in d.items():
        merged_list = value['unoccupied'] + value['occupied']
        sorted_list = sorted(merged_list, key=lambda x: x[1])
        merged_data[key] = sorted_list
    return merged_data


def get_comfort_levels(
    results: dict[str: float],
    comfort_levels: dict[str: float],
) -> dict[str: str]:
    """
    Returns a dictionnary where the keys are the traffic classes and the values are the comfort levels.
    Accelerations in 'results' are compared with the limits defined for each comfort level in the dict 'comfort_levels'.
    """
    comfort_results = {}
    for t_class, acceleration in results.items():
        for level, limit in list(comfort_levels.items()):
            if acceleration <= limit:
                comfort_results[t_class] = level
                break
    return comfort_results


# functions specific to CSA S7 procedure :

def critical_vertical_frequencies(
    bending_frequencies: list[tuple[int, float]],
    torsion_frequencies: list[tuple[int, float]],
) -> list[tuple[int, float]]:
    """
    Returns a list containing the vertical frequencies of the bridge that are in the critical range 
    as per CSA S7 Cl.7.3.2.3.
    """
    vertical_freqs = bending_frequencies + torsion_frequencies
    # sort by acending frequencies
    vertical_freqs.sort(key=lambda x: x[1])
    
    critical_freqs = [freq for freq in vertical_freqs if CRITICAL_VERTICAL_RANGE_START <= freq[1] <= CRITICAL_VERTICAL_RANGE_END]
    return critical_freqs


def critical_lateral_frequencies(
    lateral_frequencies: list[tuple[int, float]],
) -> list[tuple[int, float]]:
    """
    Returns a list containing the lateral frequencies of the bridge that are in the critical range 
    as per CSA S7 Cl.7.3.2.3.
    """     
    critical_freqs = [freq for freq in lateral_frequencies if CRITICAL_LATERAL_RANGE_START <= freq[1] <= CRITICAL_LATERAL_RANGE_END]
    return critical_freqs


def expected_vertical_walking_frequency(
    d: float,
) -> float:
    """
    Returns the expected walking frequency in the vertical direction
    as per CSA S7 Cl.7.3.3.3.
    """
    return 0.099 * d**2 - 0.644 * d + 2.188


def expected_lateral_walking_frequency(
    d: float,
) -> float:
    """
    Returns the expected walking frequency in the lateral direction
    as per CSA S7 Cl.7.3.3.3.
    """
    return  expected_vertical_walking_frequency(d) / 2


def vertical_resonating_harmonic(
    critical_freq: float,
    expected_freq: float,
) -> int:
    """
    Returns an integer corresponding to the resonating harmonic of the bridge for the vertical direction
    as per CSA S7 Cl.7.3.3.3.
    """
    return round(critical_freq / expected_freq)


def lateral_resonating_harmonic(
    critical_freq: float,
    expected_freq: float,
) -> int:
    """
    Returns an integer corresponding to the resonating harmonic of the bridge for the lateral direction
    as per CSA S7 Cl.7.3.3.3.
    """
    return round_to_nearest_odd(critical_freq / expected_freq)


def number_of_pedestrians(
    d: float,
    L: float,
    B: float,
) -> float:
    """
    Returns the number of pedestrians on the loaded surface of a pedestrian bridge.
    """
    return round(d * L * B)


def mass_pedestrians(
    d: float,
    L: float,
    B: float,
) -> float:
    """
    Returns the mass per unit length of the pedestrians on the bridge.
    
    Assumes the mass of one pedestrian is 70 kg.
    """   
    return number_of_pedestrians(d, L, B) * G0_mass / L 


def equivalent_number_of_pedestrians(
    d: float,
    L: float,
    B: float,
    damping: float
) -> float:
    """
    Returns the equivalent number of pedestrians on the bridge
    as per CSA S7 Cl.7.3.3.5.
    """
    if d < 1.0 : 
        n_eq = 10.8 * math.sqrt(damping * number_of_pedestrians(d, L, B)) / (L * B)
    else :
        n_eq = 1.85 * math.sqrt(number_of_pedestrians(d, L, B)) / (L * B)
    return n_eq


def vertical_harmonic_load(
    vertical_resonating_harmonic: int,
    equivalent_number_of_pedestrians: float,
) -> float:
    """
    Returns the uniformly distributed harmonic load amplitude in the vertical direction 
    to be considered as per CSA S7 Cl.7.3.3.5.
    Load is an area load expressed in newtons per square meters of loaded surface of bridge.
    """  
    if vertical_resonating_harmonic <= 5 :
        return G0 * VERTICAL_DYNAMIC_LOAD_FACTOR[vertical_resonating_harmonic-1] * equivalent_number_of_pedestrians
    else:
        return 0    # to avoid error in acceleration calculation


def lateral_harmonic_load(
    lateral_resonating_harmonic: int,
    equivalent_number_of_pedestrians: float,
) -> float:
    """
    Returns the uniformly distributed harmonic load amplitude in the lateral direction 
    to be considered as per CSA S7 Cl.7.3.3.5.
    Load is an area load expressed in newtons per square meters of loaded surface of bridge.
    """
    if lateral_resonating_harmonic <= 5 :
        return G0 * LATERAL_DYNAMIC_LOAD_FACTOR[lateral_resonating_harmonic-1] * equivalent_number_of_pedestrians
    else:
        return 0    # to avoid error in acceleration calculation


def S7_vertical_results(
    bridge: PedestrianBridge,
) -> dict[str: list[ tuple[int, float, float]]]:
    """
    Returns a dictionnary containing the vertical accelerations of the pedestrian bridge as per CSA S7 Cl.7.3.3.7 
    for each critical frequencies and for each traffic class defined in CSA S7 Table 7.1 for the two conditions :
        - unoccupied
        - occupied
    Assumes the bridge is a simply-supported beam with a distributed mass and a distributed harmonic loading along the length.
    """
    # store initial mass as unoccupied
    mass_unoccupied = bridge.mass
    
    vert_acc_unocc = {}
    vert_acc_occ = {}
    for (t_class, d) in S7_PED_DENSITIES.items():
        # reset mass to unoccupied
        bridge.mass = mass_unoccupied
        # calculate frequencies and accelerations for unoccupied structure
        freqs_unocc = critical_vertical_frequencies(bridge.bending_frequencies(), bridge.torsion_frequencies())
        accelerations_unocc = [(n, f_n, bridge.calculate_max_vertical_acceleration(d, bridge.mass, n, f_n)) for (n, f_n) in freqs_unocc]
        vert_acc_unocc.update({t_class: accelerations_unocc})

        # modify initial mass considering pedestrian density
        bridge.mass = mass_unoccupied + mass_pedestrians(d, bridge.L, bridge.B)
        # calculate frequencies and accelerations for occupied structure
        freqs_occ = critical_vertical_frequencies(bridge.bending_frequencies(), bridge.torsion_frequencies())
        accelerations_occ = [(n, f_n, bridge.calculate_max_vertical_acceleration(d, bridge.mass, n, f_n)) for (n, f_n) in freqs_occ]
        vert_acc_occ.update({t_class: accelerations_occ})

    # reset mass to intial mass
    bridge.mass = mass_unoccupied

    # merge both dicts into one 
    dict_to_reorder = {'unoccupied': vert_acc_unocc, 'occupied': vert_acc_occ}
    # reorder it to have 't_class' as the first level key
    dict_reordered = reorder_dict(dict_to_reorder)
    # merge results for unoccupied and occupied :
    dict_results = merge_data(dict_reordered)

    return dict_results


def S7_vertical_max_acceleration(
    bridge: PedestrianBridge,
) -> dict[str: float]:
    """
    Returns a dictionnary containing the maximum vertical acceleration of the pedestrian bridge as per CSA S7 Cl.7.3.3.7 
    for each traffic class defined in CSA S7 Table 7.1.
    Assumes the bridge is either unoccupied or occupied when calculating frequencies and maximum accelerations.
    Assumes the bridge is a simply-supported beam with a distributed mass and a distributed harmonic loading along the length.
    """
    results = S7_vertical_results(bridge) 
    
    max_acceleration = {}
    for t_class, values in results.items():
        max_value = max(values, key=lambda x: x[2])[2]
        max_acceleration[t_class] = max_value
    
    return max_acceleration


def S7_vertical_comfort(
    bridge: PedestrianBridge,
) -> dict[str: str]:
    """
    Returns a dictionnary containing the maximum vertical acceleration of the pedestrian bridge as per CSA S7 Cl.7.3.3.7 
    for each traffic class defined in CSA S7 Table 7.1.
    Assumes the bridge is either unoccupied or occupied when calculating frequencies and maximum accelerations.
    Assumes the bridge is a simply-supported beam with a distributed mass and a distributed harmonic loading along the length.
    """
    max_accelerations = S7_vertical_max_acceleration(bridge) 
    
    return get_comfort_levels(max_accelerations, COMFORT_LEVELS_VERTICAL)


def S7_lateral_results(
    bridge: PedestrianBridge,
) -> dict[str: list[ tuple[int, float, float]]]:
    """
    Returns a dictionnary containing the lateral accelerations of the pedestrian bridge as per CSA S7 Cl.7.3.3.7 
    for each critical frequencies and for each traffic class defined in CSA S7 Table 7.1 for the two conditions :
        - unoccupied
        - occupied
    Assumes the bridge is a simply-supported beam with a distributed mass and a distributed harmonic loading along the length.
    """
    # store initial mass as unoccupied
    mass_unoccupied = bridge.mass
    
    vert_acc_unocc = {}
    vert_acc_occ = {}
    for (t_class, d) in S7_PED_DENSITIES.items():
        # reset mass to unoccupied
        bridge.mass = mass_unoccupied
        # calculate frequencies and accelerations for unoccupied structure
        freqs_unocc = critical_lateral_frequencies(bridge.lateral_frequencies())
        accelerations_unocc = [(n, f_n, bridge.calculate_max_lateral_acceleration(d, bridge.mass, n, f_n)) for (n, f_n) in freqs_unocc]
        vert_acc_unocc.update({t_class: accelerations_unocc})

        # modify initial mass considering pedestrian density
        bridge.mass = mass_unoccupied + mass_pedestrians(d, bridge.L, bridge.B)
        # calculate frequencies and accelerations for occupied structure
        freqs_occ = critical_lateral_frequencies(bridge.lateral_frequencies())
        accelerations_occ = [(n, f_n, bridge.calculate_max_lateral_acceleration(d, bridge.mass, n, f_n)) for (n, f_n) in freqs_occ]
        vert_acc_occ.update({t_class: accelerations_occ})

    # reset mass to intial mass
    bridge.mass = mass_unoccupied

    # merge both dicts into one 
    dict_to_reorder = {'unoccupied': vert_acc_unocc, 'occupied': vert_acc_occ}
    # reorder it to have 't_class' as the first level key
    dict_reordered = reorder_dict(dict_to_reorder)
    # merge results for unoccupied and occupied :
    dict_results = merge_data(dict_reordered)

    return dict_results


def S7_lateral_max_acceleration(
    bridge: PedestrianBridge,
) -> dict[str: float]:
    """
    Returns a dictionnary containing the maximum lateral acceleration of the pedestrian bridge as per CSA S7 Cl.7.3.3.7 
    for each traffic class defined in CSA S7 Table 7.1.
    Assumes the bridge is either unoccupied or occupied when calculating frequencies and maximum accelerations.
    Assumes the bridge is a simply-supported beam with a distributed mass and a distributed harmonic loading along the length.
    """
    results = S7_lateral_results(bridge) 
    
    max_acceleration = {}
    for t_class, values in results.items():
        max_value = max(values, key=lambda x: x[2])[2]
        max_acceleration[t_class] = max_value
    
    return max_acceleration


def S7_lateral_comfort(
    bridge: PedestrianBridge,
) -> dict[str: str]:
    """
    Returns a dictionnary containing the maximum lateral acceleration of the pedestrian bridge as per CSA S7 Cl.7.3.3.7 
    for each traffic class defined in CSA S7 Table 7.1.
    Assumes the bridge is either unoccupied or occupied when calculating frequencies and maximum accelerations.
    Assumes the bridge is a simply-supported beam with a distributed mass and a distributed harmonic loading along the length.
    """
    max_accelerations = S7_lateral_max_acceleration(bridge) 
    
    return get_comfort_levels(max_accelerations, COMFORT_LEVELS_LATERAL)