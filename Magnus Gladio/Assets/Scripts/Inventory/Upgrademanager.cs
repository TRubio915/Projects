using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using TMPro;

public class Upgrademanager : MonoBehaviour
{
    private int dashLvl;
    public TMP_Text dashText;
    private int dashGem;
    private int mvSpeedLvl;
    public TMP_Text mvSpeedText;
    private int mvSpeedGem;
    private int dmgLvl;
    public TMP_Text dmgText;
    private int dmgGem;
    public GameObject gemCount;
    public GameObject Hint;
    public TMP_Text hint;

    public float timeToHideHint;
    private float hideHintCD;

    private GameManager manager;

    public GameObject Player;


    public GameObject MoveSpeedUpgradeButton, DamageUpgradeButton;

    void Start()
    {
        dashLvl = 0;
        mvSpeedLvl = 0;
        dmgLvl = 0;
        mvSpeedGem = 10;
        dashGem = 10;
        dmgGem = 10;

        manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
    }

    // Update is called once per frame
    void Update()
    {
        if (this.enabled)
        {
            if (hideHintCD > 0)
            {
                hideHintCD -= Time.unscaledDeltaTime;
            }

            if (hideHintCD <= 0 && Hint.activeInHierarchy)
            {
                Hint.SetActive(false);
            }

            UpgradeDescription();

            if (dmgLvl == 3)
            {
                DamageUpgradeButton.SetActive(false);
            }
            if (mvSpeedLvl == 3)
            {
                MoveSpeedUpgradeButton.SetActive(false);
            }
        }
    }

    private void UpgradeDescription()
    {
        if (dashLvl == 0)
        {
            dashText.text = "";
        }
        else if (dashLvl == 1)
        {
            dashText.text = "";
        }
        else if (dashLvl == 2)
        {
            dashText.text = "";
        }
        else if (dashLvl == 3)
        {
            dashText.text = "";
        }

        if (dmgLvl == 0)
        {
            dmgText.text = "Current Level: 0 \n\n+ 0 damage\n\nNext Level: 1\n\n+ 10 damage\n\nRequires " + dmgGem + " Gems";
        }
        else if (dmgLvl == 1)
        {
            dmgText.text = "Current Level: 1 \n\n+ 10 damage\n\nNext Level: 2\n\n+ 20 damage\n\nRequires " + dmgGem + " Gems";
        }
        else if (dmgLvl == 2)
        {
            dmgText.text = "Current Level: 2 \n\n+ 20 damage\n\nNext Level: 3\n\n+ 30 damage\n\nRequires " + dmgGem + " Gems";
        }
        else if (dmgLvl == 3)
        {
            dmgText.text = "Current Level: 3  \n\n+ 30 DMG \n\n Max Level";
        }

        if (mvSpeedLvl == 0)
        {
            mvSpeedText.text = "Current Level: 0 \n\n+ 0 Speed\n\nNext Level: 1\n\n+ 10 Speed\n\nRequires " + mvSpeedGem + " Gems";
        }
        else if (mvSpeedLvl == 1)
        {
            mvSpeedText.text = "Current Level: 1 \n\n+ 10 Speed\n\nNext Level: 2\n\n+ 20 Speed\n\nRequires " + mvSpeedGem + " Gems";
        }
        else if (mvSpeedLvl == 2)
        {
            mvSpeedText.text = "Current Level: 2 \n\n+ 20 Speed\n\nNext Level: 3\n\n+ 30 Speed\n\nRequires " + mvSpeedGem + " Gems";
        }
        else if (mvSpeedLvl == 3)
        {
            mvSpeedText.text = "Current Level: 3 \n\n+ 30 Speed \n\n Max Level";
        }
    }

    public void dashUpgrade()
    {
        
        if (dashLvl <= 2)
        {
           if (dashLvl == 1)
            {
                dashGem = 30;
            }
            else if (dashLvl == 2)
            {
                dashGem = 50;
            }
            if (manager.coinsCollected >= dashGem)
            {
                dashLvl += 1;
                manager.coinsCollected -= dashGem;
                hint.text = "Upgraded!";
            }
            else
            {
                hint.text = "Not enough Gems!";
            }
            Hint.SetActive(true);
        }

        hideHintCD = timeToHideHint;
    }
    public void mvSpeedUpgrade()
    {
        if (mvSpeedLvl <= 2)
        {
            if (manager.coinsCollected >= mvSpeedGem)
            {
                mvSpeedLvl += 1;
                manager.UpdateCoin(manager.coinsCollected - mvSpeedGem);
                Player.GetComponent<PlayerController>().speedUpgrade();
                hint.text = "Upgraded!";
                if (mvSpeedLvl == 1)
                {
                    mvSpeedGem = 30;
                }               
                else if (mvSpeedLvl == 2)
                {
                    mvSpeedGem = 50;
                }

            }
            else
            {
                hint.text = "Not enough Gems!";
            }
            
            Hint.SetActive(true);
        }
        hideHintCD = timeToHideHint;
    }
    public void dmgUpgrade()
    {
        if (dmgLvl <= 2)
        {
            if (manager.coinsCollected >= dmgGem)
            {
                dmgLvl += 1;
                Player.GetComponent<MeleeAttackController>().dmgUpgrade();
                manager.UpdateCoin(manager.coinsCollected - dmgGem);
                hint.text = "Upgraded!";

                if (dmgLvl == 1)
                {
                    dmgGem = 30;
                }
                else if (dmgLvl == 2)
                {
                    dmgGem = 50;
                }
            }
            else
            {
                hint.text = "Not enough Gems!";
            }
          
            Hint.SetActive(true);
        }
        hideHintCD = timeToHideHint;

    }

    public void HideHint()
    {
        Hint.SetActive(false);
    }

}
