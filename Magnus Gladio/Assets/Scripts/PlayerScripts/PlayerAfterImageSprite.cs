using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PlayerAfterImageSprite : MonoBehaviour
{
    private Transform player;
    private GameObject playerSpriteObject;
    private SpriteRenderer sr;
    private SpriteRenderer playerSR;
    private Color color;
    [SerializeField]
    private float activeTime = 0.1f;
    private float timeActivated;
    private float alpha;
    [SerializeField]
    private float alphaSet = 0.8f;
    [SerializeField]
    private float alphaDecay;
    private bool isPlayerDetected = false; 

    private void OnEnable()
    {
       

    }
    

    private void Update()
    {
        while(!isPlayerDetected)
        {
            sr = GetComponent<SpriteRenderer>();
            player = GameObject.FindGameObjectWithTag("Player").GetComponent<Transform>();
            playerSpriteObject = GameObject.FindGameObjectWithTag("Player");
            playerSR = playerSpriteObject.GetComponent<SpriteRenderer>();

            alpha = alphaSet;

            sr.sprite = playerSR.sprite;

            transform.position = player.position;

            transform.rotation = player.rotation;
            timeActivated = Time.deltaTime/Time.timeScale;
            isPlayerDetected = true;
        }

        alpha -= alphaDecay * Time.deltaTime / Time.timeScale;
        color = new Color(1, 1, 1, alpha);
        sr.color = color;

        if(Time.time >= timeActivated + activeTime)
        {
            //PlayerAfterImagePool.Instance.AddToPool(gameObject);
            //Add to pool
        }
    }
}
