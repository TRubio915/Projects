using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class PlayerStats
{

    public static float PlayerHealth;
    public static float PlayerStamina;
    public static float PlayerMaxHealth;
    public static float PlayerMaxStamina;


    public static void UpdatePlayerHealth(float value)
    {
        PlayerHealth = value;
    }
    public static void UpdatePlayerStamina(float value)
    {
        PlayerStamina = value;
    }

    public static void UpdateMaxPlayerHealth(float value)
    {
        PlayerMaxHealth = value;
    }
    public static void UpdateMaxPlayerStamina(float value)
    {
        PlayerMaxStamina = value;
    }


}
