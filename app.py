import streamlit as st
import plotly.graph_objects as go
import math
import app_S7_vibrations.calcs as calcs

st.title("Design of Pedestrian-Induced Vibrations according to CSA S7")

st.subheader('Assumptions:')
st.markdown('- The bridge is simply-supported.\n- The dynamic response is similar to a beam with a uniformly distributed mass and subjected to a uniformly distributed harmonic loading.')

st.sidebar.header("Bridge Parameters")
L_value = st.sidebar.number_input("Span ($m$)", value=38.85, min_value = 1.0, step= 1.0)
B_value = st.sidebar.number_input("Width accessible to pedestrians ($m$)", value=2.5, min_value = 0.5, step= 0.5)
damping_value = st.sidebar.number_input("Structural damping (%)", value=0.6, min_value = 0.05, step= 0.05 )
mass_value = st.sidebar.number_input("Mass ($kg/m$)", value=1456, min_value = 100, step= 100)
E_value = st.sidebar.number_input("Modulus of elasticity of material ($MPa$)", value=210000, min_value = 1000, step= 10000)
I_vert_value = st.sidebar.number_input("Second moment of area in the vertical direction ($m^4$)", value=0.03, min_value = 0.0001, step= 0.005)
I_lat_value = st.sidebar.number_input("Second moment of area in the lateral direction ($m^4$)", value=0.035, min_value = 0.0001, step= 0.005)
# G_value = st.sidebar.number_input("Shear modulus of material ($MPa$)", value=10000)
# J_value = st.sidebar.number_input("Second polar moment of area ($m^4$)", value=9000)
# mmi_value = st.sidebar.number_input("Mass moment of inertia ($kg.m^2$)", value=2000)


bridge = calcs.PedestrianBridge(
    L = L_value,      
    B = B_value,             
    damping = damping_value/100,        
    mass = mass_value,       
    E = E_value * 10**6,   
    I_vert = I_vert_value, 
    I_lat = I_lat_value,            
    # G = G_value * 10**6,              
    # J = J_value,              
    # mmi = mmi_value,
)

vert_results = calcs.S7_vertical_results(bridge)
vert_comfort = calcs.S7_vertical_comfort(bridge)

# create figure for vertical accelerations
fig1 = go.Figure()


# Add acceleration results
for t_class in vert_results:
    x_values = [result[1] for result in vert_results[t_class]]
    y_values = [result[2] for result in vert_results[t_class]]
    # Plot lines
    fig1.add_trace(
        go.Scatter(
        x=x_values, 
        y=y_values,
        mode = 'markers',
        marker=dict(size=10, 
                    line=dict(width=2,
                              color='DarkSlateGrey')),
        name=f"{t_class} ({calcs.S7_PED_DENSITIES[t_class]} ped/m^2)",
        )
    )

COMFORT_COLORS = ["LightSkyBlue", "LightSeagreen", "Lemonchiffon", "Firebrick"]
# Set lower limit to 0 for the first color
previous_limit = 0

# Add shapes and legends for different comfort levels
for idx, (comfort, limit) in enumerate(calcs.COMFORT_LEVELS_VERTICAL.items()):

    # Set a large value if limit is infinity for the last comfort level
    if limit == float('inf'):
        limit = previous_limit + max(y_values)

    fig1.add_shape(
        type="rect",
        xref="x",
        yref="y",
        x0=1,
        x1=13,
        y0=previous_limit,
        y1=limit,
        fillcolor=COMFORT_COLORS[idx],
        opacity=0.3,
        line=dict(width=0),
    )
    # Add invisible scatter plot as legend for the different comfort levels
    fig1.add_trace(go.Scatter(x=[None], 
                             y=[None],
                             mode='markers',
                             marker=dict(symbol="square",
                                         size=20,
                                         color=COMFORT_COLORS[idx]),
                             showlegend=True,
                             opacity=0.3,
                             name=f'{comfort} Comfort'))
    previous_limit = limit

# Add vertical lines for critical frequency limits
fig1.add_shape(
    type="line",
    x0=1,
    x1=1,
    y0=0,
    y1=limit,
    line=dict(
        color='DarkSlateGrey',
        width=3,
    )
)
fig1.add_shape(
    type="line",
    x0=13,
    x1=13,
    y0=0,
    y1=limit,
    line=dict(
        color='DarkSlateGrey',
        width=3,
    )
)


# Set minimal limits for x and y axis
fig1.update_xaxes(type="log", range=[0,math.log10(13)], )
fig1.update_yaxes(range=[0,max(y_values)+1])

fig1.layout.title.text = "Bridge Vertical Maximum Accelerations"
fig1.layout.xaxis.title = "Critical Vertical Frequency (Hz)"
fig1.layout.yaxis.title = "Acceleration (m/s^2)"

col1, col2 = st.columns(2, gap="medium")

st.plotly_chart(fig1)

for t_class, comfort in vert_comfort.items():
    st.markdown(f"{t_class} : {comfort} Comfort")