using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Spell : MonoBehaviour
{

    PlayerController playerController;

    public float speed;
    public float damage;

    private void Awake()
    {
        playerController = GameObject.FindGameObjectWithTag("Player").GetComponent<PlayerController>();

        if (playerController.GetMovementDirectionOfPlayer() < 1)
            GetComponent<SpriteRenderer>().flipX = true;

        GetComponent<Rigidbody2D>().velocity = new Vector2(playerController.GetMovementDirectionOfPlayer() * speed, 0);
    }
    private void Update()
    {
        if (GetComponent<Renderer>().isVisible == false)
        {
            DeleteObject();
        }
    }

    private void DisableCollider()
    {
        GetComponent<CircleCollider2D>().enabled = false;
    }

    private void DeleteObject()
    {
        Destroy(gameObject);
    }


    private void OnTriggerEnter2D(Collider2D collision)
    {
        
        if (collision.CompareTag("Enemy"))
        {
            float[] details = new float[2];
            details[0] = damage;
            details[1] = transform.position.x;            
            collision.SendMessage("Damage", details);
            //Destroy(gameObject);
        }
        else
        {
            Destroy(gameObject);
        }
        
    }

}
