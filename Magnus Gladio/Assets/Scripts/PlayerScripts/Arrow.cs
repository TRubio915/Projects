using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Arrow : MonoBehaviour
{

    Rigidbody2D rb;
    public bool hasHit = false;

    private Vector2 lastVel;
    public bool isCalledBack = false;

    [SerializeField]
    private float damage;

    [SerializeField]
    private LayerMask whatIsSolid;

    private RangedAttackWithTeleportation r_controller;
    // Start is called before the first frame update
    void Start()
    {
        rb = GetComponent<Rigidbody2D>();
        r_controller = GameObject.FindGameObjectWithTag("Player").GetComponent<RangedAttackWithTeleportation>();
    }

    // Update is called once per frame
    void Update()
    {
       
        if (r_controller.currentMode == 1 && hasHit == false)
        {
            float angle = Mathf.Atan2(rb.velocity.y, rb.velocity.x) * Mathf.Rad2Deg - 90;
            transform.rotation = Quaternion.AngleAxis(angle, Vector3.forward);
        }

        if(GetComponent<Renderer>().isVisible == false)
        {
            GetComponent<Rigidbody2D>().velocity = Vector2.zero;
            r_controller.resetThrowingVars();            

        }
        lastVel = GetComponent<Rigidbody2D>().velocity;
    }


    private void OnCollisionEnter2D(Collision2D collision)
    {
        
        if (collision.gameObject.tag != "Player") {
            
            isCalledBack = true;            

            if (collision.gameObject.tag != "Enemy" && collision.gameObject.tag != "Bounds")
            {
                if (collision.gameObject.tag == "Metal")
                {

                    Vector2 inDirection = lastVel;
                    float speed = inDirection.magnitude;
                    Vector2 inNormal = collision.contacts[0].normal;
                    Vector2 newVelocity = Vector2.Reflect(inDirection, inNormal);

                    GetComponent<Rigidbody2D>().velocity = newVelocity;

                }
                else
                {
                    rb.velocity = Vector2.zero;
                    rb.isKinematic = true;
                    GetComponent<Collider2D>().enabled = false;
                    hasHit = true;
                }
            }
            else
            {
                Time.timeScale = 1.0f;
                float[] attackDetails = new float[2];
                attackDetails[0] = damage;
                attackDetails[1] = transform.position.x;
                //collision.gameObject.GetComponent<Enemy>().Damage(attackDetails);
                collision.gameObject.GetComponent<Enemy>().Stun();
                Destroy(gameObject, 0.1f);
                r_controller.resetThrowingVars();                               
                //Damage Enemy
            }
        }
    }

}
