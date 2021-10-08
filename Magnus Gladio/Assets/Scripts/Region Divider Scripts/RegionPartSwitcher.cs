using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using UnityEngine;

public class RegionPartSwitcher : MonoBehaviour
{

    GameManager manager;

    public int[] regionParts;
    private bool cameraUpdated = false;

    private PlayerController playerController;

    // Start is called before the first frame update
    void Start()
    {
        manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
        playerController = GameObject.FindGameObjectWithTag("Player").GetComponent<PlayerController>();
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    private void OnTriggerEnter2D(Collider2D collider)
    {
        if (collider.gameObject.CompareTag("Player"))
        {
            int currentPart = manager.GetRegionPart();


            if (cameraUpdated == false)
            {
                int nextPart = currentPart;
                if (currentPart == regionParts[0])
                {
                    nextPart = regionParts[1];
                }
                else
                {
                    nextPart = regionParts[0];
                }
                cameraUpdated = true;
                playerController.StopAllBasicMovement();
                playerController.GetComponent<Rigidbody2D>().velocity = Vector2.zero;
                playerController.GetComponent<RangedAttackWithTeleportation>().canThrow = false;
                manager.switchCameraBoundsBetweenRegion(nextPart);
               
                GameObject.FindGameObjectWithTag("Player").GetComponent<Rigidbody2D>().AddForce(playerController.transform.right * 3000);

                playerController.SetAllBasicMovement();
                playerController.GetComponent<RangedAttackWithTeleportation>().canThrow = true;                
                
            }
        }
        if(collider.gameObject.CompareTag("Projectile"))
        {
            int currentPart = manager.GetRegionPart();
            
            int nextPart = currentPart;
            if (currentPart == regionParts[0])
            {
                nextPart = regionParts[1];
            }
            else
            {
                nextPart = regionParts[0];
            }

            playerController.GetComponent<RangedAttackWithTeleportation>().TeleportAndSwitchCamera(nextPart);
            
        }
    }
    

    private void OnTriggerExit2D(Collider2D other)
    {
        if (other.gameObject.CompareTag("Player"))
        {           
            cameraUpdated = false;           
        }
    }
}
