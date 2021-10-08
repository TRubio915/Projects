using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class EnemyAfterImageSprite : MonoBehaviour
{
    private Transform enemy;
    private GameObject enemySpriteObject;
    private SpriteRenderer sr;
    private SpriteRenderer enemySR;
    private Color color;
    [SerializeField]
    private float activeTime = 0.1f;
    private float timeActivated;
    private float alpha;
    [SerializeField]
    private float alphaSet = 0.8f;
    [SerializeField]
    private float alphaDecay;
    private bool isEnemyDetected = false; 

    private void OnEnable()
    {
       

    }
    

    private void Update()
    {
        while(!isEnemyDetected)
        {
            sr = GetComponent<SpriteRenderer>();
            enemy = GameObject.FindObjectOfType<WizardBoss>().GetComponent<Transform>();
            enemySpriteObject = GameObject.FindObjectOfType<WizardBoss>().gameObject;
            enemySR = enemySpriteObject.GetComponent<SpriteRenderer>();

            alpha = alphaSet;

            sr.sprite = enemySR.sprite;

            transform.position = enemy.position;

            transform.rotation = enemy.rotation;
            timeActivated = Time.deltaTime/Time.timeScale;
            isEnemyDetected = true;
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
