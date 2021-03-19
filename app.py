from flask import Flask
from datetime import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import *
import numpy as np
from numpy import *
from sympy import *
from sympy import symbols
#import scipy.integrate as odeint
#from scipy.integrate import quad
#from py_expression_eval import Parser
import math
app = Flask(__name__)

@app.route('/')
def homepage():
    the_time = datetime.now().strftime("%A, %d %b %Y %l:%M %p")

    return ""


    # inputed panel areas list
    panel_num_list = [1, 3, 5]

    # this is the number of seconds that we are evaluating at for delta T (12 hours)
    evaluate_num = 3600*12 # seconds in hour times number of hours

    # resulting panel area list - each of these will be the "x's" in the function
    single_panel_area = 2
    panel_area_list =[]
    for item in range(0, len(panel_num_list)):
        new_panel_area = panel_num_list[item] * single_panel_area
        panel_area_list.append(new_panel_area)

    print("printing panel_area_list")
    print(panel_area_list)
    print("")

    # # initial variable values
    slope = -25.1 # slope
    T_air = 303 # air temperature, units: K
    y_intercept = 4100 # y-intercept
    k = 0.022 # conductivity of polyurethane, units: W/(m*K)
    surface_area = 6 # of water tank, units: m^2
    thickness = 0.03 # depth/thickness of polyurethane, untis: m
    initial_water_temp = 300 # initial water temp, units: K

    # define constants for mass and heat capacity
    mass = 1000 #kg for 10k
    heat_capacity = 4179.6 # check units #

    # THE X OF THE FUNCTION IS CURRENT WATER TEMP, WHICH THEN BECOMES TIME
    x = symbols('x')

    # the x is time, y is temperature
    def glazed_plate_function(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area) -> float:
        a = (slope*panel_area) - ((k*surface_area)/thickness)
        b = (slope*-1*T_air*panel_area) + (y_intercept*panel_area) + (((k*surface_area)/thickness)*T_air)
        constant = 298*a + b
        time_expression = (1/a)*(constant*e**((x*a)/(mass*heat_capacity)) - b)
        return time_expression

    def cooling_function(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area) -> float:
        glazed_plate_evaluation = glazed_plate_function(evaluate_num, T_air, k, surface_area, thickness, initial_water_temp, panel_area)
        constant = glazed_plate_evaluation - T_air
        print(constant)
        a = ((k*surface_area)/thickness)/(mass*heat_capacity)
        print(a)
        time_expression = (constant * (e**(-1*((x-evaluate_num)*a)))) + T_air
        print(time_expression)
        return time_expression


    coil_efficiency = 1
    max_power = 1000

    def glazed_piecewise(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area):    
        glazed_time_expression = Piecewise((glazed_plate_function(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area), x <= evaluate_num), (cooling_function(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area), x > evaluate_num))
        return glazed_time_expression


    # def PV_panel_function(x, surface_area, thickness, solar_output, max_power):
    #     a = (k*surface_area)/thickness
    #     b = max_power * solar_output * coil_efficiency
    #     constant = 10 # need to figure this out
    #     PV_time_expression = T_air + (constant * e**(-a*x) + b)/(a)

    # loop through each of the panel areas and get a resulting list of expressions of time for each plate area
    def new_time_expression_list(x, T_air, k, surface_area, thickness, initial_water_temp):
        time_expression_list = []
        for item in range(0, len(panel_num_list)):
            glazed_new_time_exp = glazed_piecewise(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area_list[item])
            time_expression_list.append(glazed_new_time_exp)
        return time_expression_list

    time_expression_output = new_time_expression_list(x, T_air, k, surface_area, thickness, initial_water_temp)


    print("")
    print("printing items in time_expression_list")
    print(time_expression_output)


    # calculate delta_T for each plate area
    def new_delta_T_list(time_expression_list):
        delta_T_list = []
        at_final_time = new_time_expression_list(evaluate_num*2, T_air, k, surface_area, thickness, initial_water_temp)
        at_initial_time = new_time_expression_list(0, T_air, k, surface_area, thickness, initial_water_temp)
        for index in range(0, len(panel_num_list)):
            delta_T_list.append(at_final_time[index] - at_initial_time[index])
        return delta_T_list

    delta_T_output = new_delta_T_list(time_expression_output)

    print("")
    print("printing delta_T_list")
    print(delta_T_output)
    print("")

    # **** TO DO LATER ****
    # define function for PV panels
    # integrate function for PV panels

    # TOTAL COST OF SYSTEM
    # define constants for cost
    fixed_costs = 3000 # need a real estimate here - these are fake numbers
    panel_cost = 400 # also need a real estimate

    # function to calculate cost of the system
    def calculate_system_cost(panel_cost, panel_num, fixed_costs):
        total_cost = panel_cost * panel_num + fixed_costs
        return total_cost

    # loop through all of the panel numbers to get a list of total costs
    def new_total_cost_list():
        total_cost_list = []
        for item in panel_num_list:
            new_cost = calculate_system_cost(panel_cost, item, fixed_costs)
            total_cost_list.append(new_cost)
        
        return total_cost_list

    total_cost_output = new_total_cost_list()

    print("printing total cost list")
    print(total_cost_output)
    print("")

    # ANNUAL SAVINGS
    total_firewood_cost = 6000  # in dollars
    fraction_wood_used = 2/3    # fraction of wood that can be adjusted
    final_water_temp = 373      # in Kelvin

    # function to calculate the annual savings
    def calculate_annual_savings(delta_T, fraction_wood_used):
        annual_savings = fraction_wood_used * total_firewood_cost * (delta_T)/final_water_temp
        return annual_savings

    # loop through all the delta Ts to get a list of annual savings for each delta T
    def new_annual_savings_list(delta_T_list):
        annual_savings_list = []
        for item in delta_T_list:
            new_annual_savings = calculate_annual_savings(item, fraction_wood_used)
            annual_savings_list.append(new_annual_savings)
        
        return annual_savings_list

    annual_savings_output = new_annual_savings_list(delta_T_output)

    print("printing annual_savings_list")
    print(annual_savings_output)
    print("")

    # PAYBACK PERIOD
    # function to calculate payback period
    def calculate_payback_period(total_cost, annual_savings):
        payback_period = total_cost/annual_savings
        return payback_period

    # loop through all the total costs to get a list of payback periods for each panel number
    def new_payback_period_list(total_cost_list, annual_savings_list):
        payback_period_list = []    
        for item in range(0, len(panel_num_list)):            
            new_payback = calculate_payback_period(total_cost_list[item], annual_savings_list[item])
            payback_period_list.append(new_payback)
        return payback_period_list

    payback_period_output = new_payback_period_list(total_cost_output, annual_savings_output)

    print("printing payback_period_list")
    print(payback_period_output)
    print("")

    print("- - - - - - - -")

    complete_array = np.array([panel_num_list, panel_area_list, delta_T_output, total_cost_output, annual_savings_output, payback_period_output])

    # create figure for plots
    fig = plt.figure(figsize=(10, 6))

    #TO DO LATER
    # the percent reduction in firewood slider determines which values are displayed on the payback graph


    def plot_time_temp_graph(T_air, k, surface_area, thickness, initial_water_temp):
        # axes for temp vs time graph
        sub1_x_axis_0 = 0
        sub1_x_axis_f = evaluate_num * 2  # this should be 12 hours, in hours
        sub1_y_axis_0 = 280
        sub1_y_axis_f = 460  # boiling temp in K

        sub1 = plt.subplot(2, 2, 1)
        sub1.cla()  # clear the subplot's axes before drawing on it
        sub1.set_title('temp vs. time')
        # sub1.set_xticks(np.arange(sub1_x_axis_0, sub1_x_axis_f, step=3600))
        # sub1.set_yticks(np.arange(sub1_y_axis_0, sub1_y_axis_f, step=30))
        sub1.set_xlabel('time (sec)')
        sub1.set_ylabel('temperature (K)')
        
        # x = np.linspace(sub1_x_axis_0, sub1_x_axis_f, 200)
        x = np.linspace(0, evaluate_num * 2, 200)
        
        sub1.set_xlim([sub1_x_axis_0, sub1_x_axis_f])
        sub1.set_ylim([sub1_y_axis_0, sub1_y_axis_f])

        temp_time_graph_list = []

        for item in range(0, len(panel_num_list)):
            pw = np.piecewise(x, [x < evaluate_num, x >= evaluate_num], [
                lambda z: glazed_plate_function(z, T_air, k, surface_area, thickness, initial_water_temp,
                                                panel_area_list[item]), lambda q: cooling_function(q, T_air, k, surface_area, thickness, initial_water_temp, panel_area_list[item])])
            temp_time_graph_list.append(sub1.plot(x, pw))




    plot_time_temp_graph(T_air, k, surface_area, thickness, initial_water_temp)

    def plot_payback_graph(delta_T_output, payback_period_output):
        # axes for payback period vs delta T graph
        sub2_x_axis_0 = 0
        sub2_x_axis_f = 373.16  # boiling temp in K
        sub2_y_axis_0 = 0
        sub2_y_axis_f = 10     # hopefully not more than 10 years

        sub2 = plt.subplot(2, 2, 2)
        sub2.set_title('payback vs. delta T')
        sub2.set_xticks(np.arange(sub2_x_axis_0, sub2_x_axis_f, step=30))
        sub2.set_yticks(np.arange(sub2_y_axis_0, sub2_y_axis_f, step=1))
        sub2.set_xlabel('payback period (yrs)')
        sub2.set_ylabel('delta T (K)')
        #sub2.set_xlim([sub2_x_axis_0, sub2_x_axis_f]) 
        #sub2.set_ylim([sub2_y_axis_0, sub2_y_axis_f])
        f_d2, = sub2.plot(payback_period_output, delta_T_output)
        # f_d2, = sub2.plot(delta_T_output, payback_period_output)

    plot_payback_graph(delta_T_output, payback_period_output)

    # SLIDERS
    # Create axes for sliders
    sub1_T_air = fig.add_axes([0.15, 0.34, 0.25, 0.02])
    sub1_T_air.spines['top'].set_visible(True)
    sub1_T_air.spines['right'].set_visible(True)

    sub1_single_panel_area = fig.add_axes([0.15, 0.3, 0.25, 0.02])
    sub1_single_panel_area.spines['top'].set_visible(True)
    sub1_single_panel_area.spines['right'].set_visible(True)

    sub1_k = fig.add_axes([0.15, 0.26, 0.25, 0.02])
    sub1_k.spines['top'].set_visible(True)
    sub1_k.spines['right'].set_visible(True)

    sub1_surface_area = fig.add_axes([0.15, 0.22, 0.25, 0.02])
    sub1_surface_area.spines['top'].set_visible(True)
    sub1_surface_area.spines['right'].set_visible(True)

    sub1_thickness = fig.add_axes([0.15, 0.18, 0.25, 0.02])
    sub1_thickness.spines['top'].set_visible(True)
    sub1_thickness.spines['right'].set_visible(True)

    sub1_initial_water_temp = fig.add_axes([0.15, 0.14, 0.25, 0.02])
    sub1_initial_water_temp.spines['top'].set_visible(True)
    sub1_initial_water_temp.spines['right'].set_visible(True)

    sub1_fixed_costs = fig.add_axes([0.15, 0.10, 0.25, 0.02])
    sub1_fixed_costs.spines['top'].set_visible(True)
    sub1_fixed_costs.spines['right'].set_visible(True)

    sub1_panel_cost = fig.add_axes([0.15, 0.06, 0.25, 0.02])
    sub1_panel_cost.spines['top'].set_visible(True)
    sub1_panel_cost.spines['right'].set_visible(True)

    # create sliders for each of those initialized vars
    s_T_air = Slider(ax=sub1_T_air, label='air temp ', valmin=0.9*T_air, valmax=1.1*T_air, valinit=T_air, valfmt=' %2.2f K', facecolor='#006400')
    s_single_panel_area = Slider(ax=sub1_single_panel_area, label='single panel area ', valmin=0.2*single_panel_area, valmax=1.8*single_panel_area, valinit=single_panel_area, valfmt=' %2.2f m^2', facecolor='#006400')
    s_k = Slider(ax=sub1_k, label='k (conductivity) ', valmin=0.1*k, valmax=1.9*k, valinit=k, valfmt=' %2.2f W/mK', facecolor='#006400')
    s_surface_area = Slider(ax=sub1_surface_area, label='surface area ', valmin=0.1*surface_area, valmax=1.9*surface_area, valinit=surface_area, valfmt=' %2.2f m^2', facecolor='#006400')
    s_thickness = Slider(ax=sub1_thickness, label='thickness ', valmin=0.1*thickness, valmax=1.9*thickness, valinit=thickness, valfmt=' %2.2f m', facecolor='#006400')
    s_initial_water_temp = Slider(ax=sub1_initial_water_temp, label='initial water temp ', valmin=0.1*initial_water_temp, valmax=1.9*initial_water_temp, valinit=initial_water_temp, valfmt=' %2.2f K', facecolor='#006400')
    s_fixed_costs = Slider(ax=sub1_fixed_costs, label='fixed costs ', valmin=0.1*fixed_costs, valmax=1.9*fixed_costs, valinit=fixed_costs, valfmt=' %2.2f $', facecolor='#006400')
    s_panel_cost = Slider(ax=sub1_panel_cost, label='panel cost ', valmin=0.1*panel_cost, valmax=1.9*panel_cost, valinit=panel_cost, valfmt=' %2.2f $', facecolor='#006400')

    # these are more complicated
    # s_panel_number = Slider(sub1=sub1_panel_number, label='panel number ', valmin=0, valmax=10, valinit=panel_number, valfmt=' %i K', facecolor='#006400')
    # s_firewood_reduction = Slider(sub1=sub1_firewood_reduction, label='percent firewood reduction ', valmin=0, valmax=10, valinit=T1_initial, valfmt=' %i K', facecolor='#006400')




    # Update values
    def update(val):
        T_air_new = s_T_air.val
        single_panel_area_new = s_single_panel_area.val
        k_new = s_k.val
        surface_area_new = s_surface_area.val
        thickness_new = s_thickness.val
        initial_water_temp_new = s_initial_water_temp.val
        fixed_costs_new = s_fixed_costs.val
        panel_cost_new = s_panel_cost.val

        time_expression_output_new = new_time_expression_list(x, T_air_new, k_new, surface_area_new, thickness_new, initial_water_temp_new)
        
        print("")
        print("printing time_expression_output_new")
        print(time_expression_output_new)


        delta_T_output_new = new_delta_T_list(time_expression_output_new)
        total_cost_output_new = new_total_cost_list()
        annual_savings_output_new = new_annual_savings_list(delta_T_output_new)
        payback_period_output_new = new_payback_period_list(total_cost_output_new, annual_savings_output_new)

        print("")
        print("printing delta_T_list")
        print(delta_T_output_new)
        print("")


        plot_time_temp_graph(T_air_new, k_new, surface_area_new, thickness_new, initial_water_temp_new)
        plot_payback_graph(delta_T_output_new, payback_period_output_new)


    # f_d2.set_data(delta_T_output, payback_period_output)
    fig.canvas.draw_idle()
    

    s_T_air.on_changed(update)
    s_single_panel_area.on_changed(update)
    s_k.on_changed(update)
    s_surface_area.on_changed(update)
    s_thickness.on_changed(update)
    s_initial_water_temp.on_changed(update)
    s_fixed_costs.on_changed(update)
    s_panel_cost.on_changed(update)


    # lol this doesn't work
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()


    # TO DO LATER
    # create a toggle for whether glazed plate and PV panel graph is showing




if __name__ == '__main__':
    app.run(debug=True, use_reloader=True)

