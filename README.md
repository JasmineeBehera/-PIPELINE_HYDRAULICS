Project Description:
 This script models a fluid transport pipeline system to determine its
 normal operating requirements and to design a critical safety device.

Part 1: Normal Operation Hydraulics
  - Calculates pressure drop from pipe friction (Major Losses) using the
     Darcy-Weisbach equation with the Haaland approximation for the
     friction factor.
   - Calculates pressure drop from fittings (Minor Losses).
  - Determines the total required pump head and discharge pressure.

 Part 2: Safety System Design (PRV Sizing)
   - Simulates a "blocked outlet" scenario, which is a critical cause of
     overpressure in a pumped system.
   - Calculates the required orifice area for a Pressure Relief Valve (PRV)
    based on API 520 standards to safely vent the fluid.
   - Selects the next standard API orifice size.

 Relevant Subjects Utilized:
  CHC202 Fluid and Particle Mechanics: Reynolds number, friction factor, head loss.
  CHC201 Chemical Process Calculations: Mass and energy balances.
  CSI101 Computer Programming: Structuring code, creating functions.
