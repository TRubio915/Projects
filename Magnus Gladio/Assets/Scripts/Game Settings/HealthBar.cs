using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;
using UnityEngine.UI;

public class HealthBar : MonoBehaviour
{

    public Image healthCircleBar;
    public Image healthExtraBar;

    public Image healthCirclceWhite;
    public Image healthBarWhite;

    public float circlePercentage = 0.30f;

    private bool healthChanging = false;
    private bool healthBarWhiteEmpty = false;


    // Update is called once per frame
    public void UpdateHealthBar()
    {

        if (healthBarWhiteEmpty)
        {
            CircleFill();
        }
        else
        {
            ExtraBarFill();
        }
    }


    void CircleFill()
    {
        float healthPercentage = PlayerStats.PlayerHealth / PlayerStats.PlayerMaxHealth;       
        float circleFill = healthPercentage / circlePercentage;
        circleFill = Mathf.Clamp(circleFill, 0, 1);
        healthCircleBar.fillAmount = circleFill;

        if(!healthChanging)
            StartCoroutine(ChangeHealthCircleWhiteOverTime(circleFill));

    }

    IEnumerator ChangeHealthCircleWhiteOverTime(float value)
    {
        healthChanging = true;
        
        while (healthCirclceWhite.fillAmount != value)
        {

            yield return null;
            healthCirclceWhite.fillAmount = Mathf.MoveTowards(healthCirclceWhite.fillAmount, value, Time.unscaledDeltaTime/2);
        }
        healthChanging = false;
    }


    void ExtraBarFill()
    {
        float circleAmount = circlePercentage * PlayerStats.PlayerMaxHealth;
        float extraHealth = PlayerStats.PlayerHealth - circleAmount;        
        float extraTotalhealth = PlayerStats.PlayerMaxHealth - circleAmount;
        float extraFill = extraHealth/extraTotalhealth;
        extraFill = Mathf.Clamp(extraFill, 0, 1);
        healthExtraBar.fillAmount = extraFill;
        if (!healthChanging)
            StartCoroutine(ChangeHealthBarWhiteOverTime(extraFill));
    }

    IEnumerator ChangeHealthBarWhiteOverTime(float value)
    {
        healthChanging = true;
        while (healthBarWhite.fillAmount != value)
        {
            yield return null;

            healthBarWhite.fillAmount = Mathf.MoveTowards(healthBarWhite.fillAmount, value, Time.unscaledDeltaTime/2);
        }
        healthChanging = false;

        if(healthBarWhite.fillAmount <= 0)
            healthBarWhiteEmpty = true;
    }


}
