# [PyGWP](https://github.com/lfaucheux/PyGWP) - A **CO2-equivalent computer** based on static or dynamic CO2-relative global warming potentials.

## Sources 
- **[IPCC WG1 Fourth Assessment Report, 2007](https://www.ipcc.ch/publications_and_data/ar4/wg1/en/ch2s2-10-2.html)**
- **[Levasseur et al, 2010](http://pubs.acs.org/doi/abs/10.1021/es9030003)**

## Installation
    git clone git://github.com/lfaucheux/PyGWP.git
    cd PyGWP
    python setup.py install

## Requirements

- **[numpy](http://www.numpy.org/)**

## Use cases

- **Scientific modelling**

## Example usage:

    >>> from PyGWP import gwp_based_co2eq
    >>> dyn_gwp20  = gwp_based_co2eq(
    ...     first_year=2020,
    ...     project_horizon=5,
    ...     GWP_horizon=20,
    ...     static=False
    ... )
    >>> ghgs_weight_per_weight_of_output_inventory_flow = {'CO2':1.,'N2O':.0,'CH4':.0}
    >>> co2eq_traj = dyn_gwp20.co2eq_yields_trajectory_computer(
    ...     ghgs_weight_per_weight_of_output_inventory_flow,
    ...     as_row_array=False
    ... )
    >>> print(co2eq_traj['as_array'])
    [[ 1.        ]
     [ 0.95764081]
     [ 0.91469171]
     [ 0.87112496]
     [ 0.82691128]]
    >>> print(co2eq_traj['as_dict'])
    {2024: 0.82691127746144444, 2020: 1.0, 2021: 0.95764080833063492, 2022: 0.91469171438570718, 2023: 0.87112496115582216}

    >>> co2eq_traj = dyn_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':1.,'CH4':.0})
    >>> print(co2eq_traj['as_array'])                                           ^^
    [[ 292.33637282  278.90543843  265.35617058  251.68752668  237.89845498]]


    >>> co2eq_traj = dyn_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':.0,'CH4':1.})
    >>> print(co2eq_traj['as_array'])                                                    ^^
    [[ 72.2209832   70.75950679  69.17102216  67.44449179  65.56791893]]

    >>> sta_gwp20  = gwp_based_co2eq(
    ...     first_year=2020,
    ...     project_horizon=5,
    ...     GWP_horizon=20,
    ...     static=True
    ... )          ^^^^                                            
    >>> co2eq_traj = sta_gwp20.co2eq_yields_trajectory_computer({'CO2':.0,'N2O':.0,'CH4':1.})
    >>> print(co2eq_traj['as_array'])
    [[ 72.2209832  72.2209832  72.2209832  72.2209832  72.2209832]]


## License
Distributed under the [MIT license](https://opensource.org/licenses/MIT)
