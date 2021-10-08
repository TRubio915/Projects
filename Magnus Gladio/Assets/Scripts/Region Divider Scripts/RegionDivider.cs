using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RegionDivider : MonoBehaviour
{

    //public GameObject[] regions;
    private GameManager manager;

    public int[] RegionNumbers;

    private bool switchCompleted = true;

    private RangedAttackWithTeleportation r_Controller;
    

    // Start is called before the first frame update
    void Start()
    {
        manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
        r_Controller = GameObject.FindGameObjectWithTag("Player").GetComponent<RangedAttackWithTeleportation>();
    }

    
    // Update is called once per frame
    private void OnTriggerEnter2D(Collider2D other)
    {
        if (other.tag == "Player" )
        {
            if (switchCompleted == true)
            {
                //manager.setCameraBounds();
                if (manager.GetSceneNumber() == RegionNumbers[0])
                {
                    manager.setCameraBounds(RegionNumbers[1]);
                    
                }
                else
                {
                    manager.setCameraBounds(RegionNumbers[0]);
                    
                }
                GameObject.Find("SceneTransitionSprite").GetComponent<SpriteRenderer>().enabled = true;
                GameObject.Find("SceneTransitionSprite").GetComponent<Animator>().enabled = true;
                GameObject.Find("SceneTransitionSprite").GetComponent<Animator>().Play("Transition", -1, 0f);
                switchCompleted = false;
            }
            r_Controller.resetThrowingVars();
        }       

        if(other.tag == "Projectile")
        {
            if (manager.GetSceneNumber() == RegionNumbers[0])
            {
                r_Controller.TeleportAndSwitchRegion(RegionNumbers[1]);
            }
            else
            {
                r_Controller.TeleportAndSwitchRegion(RegionNumbers[0]);
            }
        }
    }
    private void OnTriggerStay2D(Collider2D other)
    {

        if (other.tag == "Player")
        {
            switchCompleted = false;
        }

    }
    private void OnTriggerExit2D(Collider2D other)
    {
        if (other.tag == "Player")
        {
            if(!switchCompleted)
            {
                switchCompleted = true;

            }
            
        }
    }
}
